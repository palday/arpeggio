#' ---
#' title: Statistical Analysis with Covariates for Musical Prediction Study
#' author: Phillip M. Alday
#' date: September 2018
#' ---

# Copyright 2018, Phillip Alday
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library("here")
library("tidyverse")  # for plotting and data manipulation
library("lme4")
library("car")
library("effects")
library("caret")
library("lattice")
library("R.matlab")
library("eegUtils") # from https://github.com/craddm/eegUtils/
library("knitr")
library("MASS") # for contr.sdiff
library("lmerOut")  # from https://bitbucket.org/palday/lmerout/
options(contrasts = c("contr.Sum","contr.Poly"))

#' columns in this data file
#'
#' - trial number
#' - subject number
#' - condition
#' - averaged voltage in the N400 window [300, 500]
#' - averaged voltage in the baseline window [-250 0]
#' - rhyme cloze probability
#' - word frequency
#' - phonological neighborhood density
#' - semantic distance from the prime word
#' - semantic distance from the lead-in sentence
#' - concreteness
#' - rhyme evaluation
#' - plausibility evaluation
#'
#' For the voltage, we have56 preprocessed channels:
#' 64 electrodes - 4 oculograms - 1 mastoid - 3 noisy occipital channels
dat <- read_csv("mlm_inputFM_prepro.csv") %>%
  # pull all the amplitude data together into long format
  gather(., key="win_chan", value="amplitude", N4_C1:BS_C60) %>%
  # split window and channel into two columns
  separate(., win_chan, c("win","chan")) %>%
  # separate out baseline and N400 windows into separate columns
  spread(., win, amplitude) %>%
  # make subject and condition into categorical variables
  mutate(subject_id = factor(subject_id),
         condition = factor(condition, levels=c(1,2,3), labels=c("congruent","intermediate","incongruent"))
         #condition = relevel(condition, "incongruent")
         )

contrasts(dat$condition)
sdiff <- contr.sdif(3)

colnames(sdiff) <- c(paste0("[",rownames(contrasts(dat$condition))[2]," > ",rownames(contrasts(dat$condition))[1],"]"),
                     paste0("[",rownames(contrasts(dat$condition))[3]," > ",rownames(contrasts(dat$condition))[2],"]"))

contrasts(dat$condition) <- sdiff
contrasts(dat$condition)

#' Load and convert channel coordinates from spherical to cartesian. For
#' formulae, see:
#'
#' - [1] https://sccn.ucsd.edu/pipermail/eeglablist/2006/001655.html
#' - [2] https://github.com/mne-tools/mne-python/blob/master/mne/transforms.py#L694-L703
#'
#' At some point, I acquired this coordinates file -- I think I got it from
#' EasyCap, the manufacturer.

channel_locations <- read.table("64Ch_actiCAP_Equidistant10_polar.txt",header=TRUE,sep="\t")

#' In Cartesian coordinates, (0,0,0) is the center of the skull and (0,0,1) is
#' the apex, i.e. electrode Cz in the 10-20 system. Negative x are left lateral,
#' positive x are right lateral. Positive y is anterior, negative y is
#' posterior.

d2r <- function(x) x * pi / 180
polar2cart <- function(theta, phi, r=1, degrees=TRUE){
  if(degrees){
    phi <- d2r(phi)
    theta <- d2r(theta)
  }
  list(x = r * sin(phi) * cos(theta),
       y = r * sin(phi) * sin(theta),
       z = r * cos(phi))
}

channel_locations <- channel_locations %>%
  as_tibble() %>%
  do(as.data.frame(polar2cart(.$theta,.$phi))) %>%
  bind_cols(channel_locations) %>%
  mutate(chan = as.character(chan))

#' For display purposes, we also need the projection into the 2D plane. I tried
#' to figure out the formulae for this and always messed it up, so we'll just
#' "steal" [the data from
#' FieldTrip](https://github.com/fieldtrip/fieldtrip/blob/master/template/layout/mpi_customized_acticap64.mat).
#'

pos <- readMat("mpi_customized_acticap64.mat")$lay[,,1]$pos
channel_plot_locations <- data.frame(x=pos[,1],y=pos[,2])[1:60,] %>%
  as_tibble() %>%
  mutate(electrode=as.character(1:60))


#' Now add in topographical information to our original data
dat <- dat %>%
  # strip out the C prefix before the channels to match our coordinates
  mutate(chan=substring(chan,2)) %>%
  left_join(channel_locations)

#' Item number is labeled `trlnumber` here. Items aren't the same as trials --
#' the same items can be presented in different orders, and it's not clear if
#' trlnumber here refers to the sequence number or the particular lexical item.
#' Items should definitely be included in the random effects the same way
#' subjects are. Although given the number of controls implemented in the
#' covariates, you should be able to get away with (1|item) instead of
#' (1+condition|item).
#'
#' We also scaled the covariates. This helps the numerical aspects and also
#' makes it easier to compare the relative weighting of the covariates.
dat <- dat %>%
  mutate_at(vars(cloze:plaus_eval), scale)

#' I'm not sure we should model both sentt_semdist and pt_semdist. They're
#' closey related conceptually and strongly correlated numerically.

cor.test(dat$sentt_semdist,dat$pt_semdist)

#' So I'll leave one out (if you have theoretical reasons to prefer one to the
#' other, feel free to change), which should also speed up model fitting!


#' I was getting some interesting warnings about rank deficiency, which often
#' happens when certain combinations don't exist / can be expressed as
#' combinations of the other columns.

X <- model.matrix(~ BS + condition * x * y * (cloze  + wordfrq + phon_nd + sentt_semdist + pt_semdist + concreteness + rhyme_eval + plaus_eval),data=dat)
lcs <- findLinearCombos(X)
colnames(X)[lcs$remove]

#' Something weird is going on with condition B and the cloze probability
#' manipulation. Is your manipulation purely cloze probability?
#'
#' Note also that there a bunch of NAs for cloze probability:

dat %>% group_by(condition) %>%
  summarize(mean=mean(cloze,na.rm=TRUE),
            sd=sd(cloze,na.rm=TRUE),
            NAs=sum(is.na(cloze)))
#' All rows containing NAs in a predictor are dropped from the model. So maybe
#' the better thing to do is omit cloze probability.
#'
#' Also, I think rhyme_eval is also an overspecification -- aren't the different
#' conditions about different rhymes?
#'
#' If so, then you could model rhyme_eval *instead of* condition. The
#' random-effects by item actually cover a lot of this type of small difference.
#' This is why there's a huge literature on how important item-level random
#' effects are in language research!

#' I've also added in the third topographical predictor (z). This is a great
#' example of issues in model fitting -- the simpler model without z didn't
#' coverge, but the more complex one did.
#'
#' I've also removed the interactions between topography and covariates: they
#' drastically increase fitting time, don't really improve fit, and greatly
#' increase model complexity, so out they go!
#'
#' This gives us the power to try testing out interactions between condition,
#' topography and the baseline interval. Adding in those interactions improves
#' model fit by all measures (AIC, BIC, LRT) so they stay.
#'
#' I've also gone ahead and scaled both EEG measures (baseline and N4) -- this
#' puts everything on the scale of standard deviations, so coefficients are
#' essentially effect sizes in Cohen's $d$. The model fit is also much better:
#' linear transformations such as scaling don't impact the model fit in terms of
#' the pure mathematics, but they can make the optimizer's job easier. With this
#' extra boost in fit, I was able to add in a by-condition slope for item.


# so that I don't have to remember the different options different optimizers have...
bobyqa <- lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e6))
nlopt <- lmerControl(optimizer="nloptwrap",optCtrl=list(maxeval=1e6))


#+ model, cache=TRUE
system.time(m <- lmer(scale(N4) ~ scale(BS) * condition * x * y * z +
                              condition * (wordfrq + phon_nd + pt_semdist + concreteness + plaus_eval) +
            (1 + condition | subject_id) +
            (1 + condition | trlnumber),
          data=dat,
          control=bobyqa,
          REML=FALSE))

#' If you don't understand some of the plots here, don't worry. The diagnostics
#' look fine overall. It looks like our residual distribution is heavy tailed
#' and more t than normal, but there's nothing to be done for that now and it
#' doesn't actually impace the inferences we really care about that much.

#+ output, cache=TRUE, results='asis'

cat(pprint(summary(m), type="html"))
# print(summary(m),correlation=FALSE,symbolic.cor=TRUE)

cat(sprintf("Number of fixed-effect correlations > 0.1: %d",sum(as.matrix(vcov(m)) > 0.1)))

#+ diagnostics, cache=TRUE, fig.width=10, fig.height=10
plot(m)

qqmath(m)

qqmath(ranef(m,condVar=TRUE))

dotplot(ranef(m,condVar=TRUE))

fortify.merMod(m, drop_na(dat)) %>%
  ggplot(aes(.fitted,.resid)) +
    geom_point(aes(color=condition),alpha=0.3) +
    facet_wrap( ~ subject_id) +
    geom_hline(yintercept=0) +
    theme_light() +
    labs(title="Residuals vs. Fitted",
         subtitle="Should be cloud shaped",
         x="Fitted",
         y="Residuals")


#+ anova, cache = TRUE
(a <- Anova(m))

#' For the ANOVA-style display, we can also omit the baseline interval because
#' we don't actually care about it, even if we have to model it.

# a[!str_detect(rownames(a),"BS"),]

#+ anova_output, results='asis'
kable(a[!str_detect(rownames(a),"BS"),])

#' The above summary can also be expressed graphically.

# Wald is the fastest method; boot the most accurate and slowest, profile is a
# the middle road
#+ cicalc, cache=TRUE

ci <- confint(m,method="Wald")

# If doing profile CIs, then it's worthwhile to use profile() explicitly and
# change some options

#system.time(ci <- profile(m,optimizer="nloptwrap",parallel="multicore",ncpus=2))

#+ coefplot, cache=TRUE, fig.width=7, fig.height=5
ci.gg <- ci %>%
  as_tibble(rownames = "coef") %>%
  dplyr::filter(substr(coef,1,4) != ".sig") %>% # omit random effects
  mutate(est = fixef(m),
         coef = factor(coef, levels=rev(names(fixef(m))))) %>%
  dplyr::filter(coef != "(Intercept)") %>% # intercept
  dplyr::filter(!str_detect(coef,"BS")) %>% # baseline terms
  mutate(topo=as.logical(str_detect(coef,"x|y|z")))


gg <- ggplot(mapping=aes(x=coef,y=est,ymin=`2.5 %`, ymax=`97.5 %`)) +
  geom_pointrange(fatten=1.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Coefficient",y="Estimate (standardized difference)") +
  coord_flip() +
  theme_light()

gg %+% subset(ci.gg, !topo) + ggtitle("Global Effects")

gg %+% subset(ci.gg, as.logical(str_detect(coef,"x"))) + ggtitle("Lateral Effects")
gg %+% subset(ci.gg, as.logical(str_detect(coef,"y"))) + ggtitle("Saggital Effects")


#' We use [`eegUtils`](https://github.com/craddm/eegUtils/) for the topoplot.
#' If you have `devtools`` installed, then install using
#' `devtools::install_github("craddm/eegUtils")`.
#'
#' First, we do the equivalent of grand-average plots in the critical time
#' window.

#+ topo-plots
dat.plot <- drop_na(dat) %>% mutate(fitted=fitted(m)) %>% as_tibble()
dat.plot <- dat.plot %>%
  group_by(chan, condition, subject_id) %>%
  summarize(amplitude=mean(N4-BS,na.rm=TRUE)) %>%
  group_by(chan, condition) %>%
  summarize(amplitude=mean(amplitude)) %>%
  ungroup() %>%
  rename(electrode=chan) %>%
  left_join(channel_plot_locations) %>%
  drop_na() %>%
  # center
  mutate(x = x-mean(x),
         y = y-mean(y)) %>%
  # scale
  mutate(x = 0.5 * x/max(abs(x)),
         y = 0.5 * y/max(abs(y)))

dat.plot %>% subset(condition == "congruent") %>% topoplot(chan_marker="name")

dat.plot %>% subset(condition == "intermediate") %>% topoplot()

dat.plot %>% subset(condition == "incongruent") %>% topoplot()

#' We can also do difference waves:

dat.plot %>% spread(condition,amplitude) %>%
  # need to specify dplyr:: because MASS also has a select
  dplyr::select(-intermediate) %>% mutate(amplitude = incongruent - congruent) %>%
  topoplot()

#' We can check that out model make the right predictions. We skip the baseline
#' correction so that we avoid some scaling issues. It's all about the
#' comparison anyway, but this will make this new plot look a bit different than
#' the previous ones.

drop_na(dat) %>% mutate(N4 = scale(N4)) %>%
  group_by(chan, condition, subject_id) %>%
  summarize(amplitude=mean(N4,na.rm=TRUE)) %>%
  group_by(chan, condition) %>%
  summarize(amplitude=mean(amplitude)) %>%
  ungroup() %>%
  rename(electrode=chan) %>%
  left_join(channel_plot_locations) %>%
  drop_na() %>%
  mutate(x = x-mean(x),
         y = y-mean(y)) %>%
  mutate(x = 0.5 * x/max(abs(x)),
         y = 0.5 * y/max(abs(y))) %>%
  subset(condition == "congruent") %>%
  topoplot(method="gam")+ ggtitle("congruent: Observed Values")

drop_na(dat) %>% mutate(fitted=fitted(m)) %>%
  group_by(chan, condition, subject_id) %>%
  summarize(amplitude=mean(fitted,na.rm=TRUE)) %>%
  group_by(chan, condition) %>%
  summarize(amplitude=mean(amplitude)) %>%
  ungroup() %>%
  rename(electrode=chan) %>%
  left_join(channel_plot_locations) %>%
  drop_na() %>%
  mutate(x = x-mean(x),
         y = y-mean(y)) %>%
  mutate(x = 0.5 * x/max(abs(x)),
         y = 0.5 * y/max(abs(y))) %>%
  subset(condition == "congruent") %>%
  topoplot(method="gam") + ggtitle("congruent: Fitted Values")

drop_na(dat) %>% mutate(residuals=residuals(m)) %>%
  group_by(chan, condition, subject_id) %>%
  summarize(amplitude=mean(residuals,na.rm=TRUE)) %>%
  group_by(chan, condition) %>%
  summarize(amplitude=mean(amplitude)) %>%
  ungroup() %>%
  rename(electrode=chan) %>%
  left_join(channel_plot_locations) %>%
  drop_na() %>%
  mutate(x = x-mean(x),
         y = y-mean(y)) %>%
  mutate(x = 0.5 * x/max(abs(x)),
         y = 0.5 * y/max(abs(y))) %>%
  subset(condition == "congruent") %>%
  topoplot(method="gam") + ggtitle("congruent: Residual Values")

#' Hmmm, the topographical fits aren't perfect, but they're okay. Remember that
#' all of this is on the standard deviation / unit scale and not the original
#' scale.

#    spread(condition,amplitude) %>%
#    select(-B) %>%
#    mutate(amplitude = A - C) %>%
#    topoplot()

#' We could achieve something similar using the `eegkit` package and this command:
#' `eegspace(dat.plot[,c("x","y")]*20,dat.plot$amplitude)`.


#' For functions that applied within the model formula, we can often pull out
#' aspects of their computation. This is useful for undoing scaling, which we will use below.
unscale_params <- function(model, term){
list(center=attr(model@frame[[sprintf("scale(%s)",term)]],"scaled:center"),
     scale=attr(model@frame[[sprintf("scale(%s)",term)]],"scaled:scale"))
}


#+ effects, cache=TRUE
system.time(e <- allEffects(m, KR=FALSE, xlevels=7))

#' We're limited in how much we can adapt the default plots without a lot of
#' effort so we just convert to a dataframe and use our old friend ggplot2. For
#' the covariates, we don't have to worry about topography, because they don't
#' interact with topographical factors.

# skip the final x-y stuff
#plot(e[1:5],rug=FALSE,multiline=TRUE,ci.style="band")

#+ effects-plot, cache=TRUE, fig.width=10, fig.height=3
# convert the list of effects to list of dataframes
edf <- lapply(e[1:5],as.data.frame)
# add in a column for each dataframe containing the name of the effect
edf <- Map(cbind, edf, effect = names(edf))
# combine and plot
# we can ignore the warning about converstion to character for factors
edf %>% bind_rows() %>% as_tibble() %>%
  # get into long format so that we can treat each covariate as a facet/panel
  gather(key="covariate",value="val",
         -effect, -condition, -fit, -se, -lower, -upper) %>%
  ggplot(aes(x=val,y=fit,ymin=lower,ymax=upper,color=condition,fill=condition)) +
  geom_line() +
  # we don't want the out edges of the ribbons marked
  geom_ribbon(color=NA,alpha=0.3) +
  # free_x because the standardized effects are still on different scales
  facet_grid( ~ effect,
              scales="free_x",
              labeller=labeller(effect=function(x) {sub("condition:","",x,fixed=TRUE)})) +
  theme_light() +
  scale_y_reverse() + # plot negativity upward (not unusual in ERP)
                      # but also "up" is now a bigger N400 effect
  labs(x="Standardized Units",
       y="Amplitude (standard deviations, inverse scale)",
       title="Effects as modelled",
       subtitle="with 95% confidence interval") -> g.eff.scaled
print(g.eff.scaled)

shift <- unscale_params(m,"N4")$center
stretch <- unscale_params(m,"N4")$scale
unscale <- function(x) (x * stretch) + shift

g.eff.scaled +
  aes(y=unscale(fit),ymin=unscale(lower),ymax=unscale(upper)) +
  labs(y="Amplitude (µV, inverse scale)")

#' Overall, this looks good. We see that the different covariates have some
#' moderating influence on the effect of condition, but none really change the
#' overall structure, as the lines are largely parallel-ish and the overall
#' (vertical) order of the effects doesn't change  (within the uncertainty given by the
#' confidence intervals)

#' Note wordfrq looks to have different scale, even on the standardized scale.
#' This because of a long right tail: you have a few very high frequency
#' words, while the majority are all about the same frequency.

plot(density(dat$wordfrq),main="Distribution of Word Frequencies")

#' (NB for this plot:  Standardization doesn't impact the overall shape of the
#' density plot -- it translates it (slides it along the x-axis) and can squish
#' or stretch out the absolute values on the x-axis, but doesn't impact the
#' relative values.)

# e.con <- as.data.frame(e[["condition:x:y:z:concreteness"]])
# e.con %>% group_by(condition,concreteness,x,y) %>%
#   summarise_all(mean) %>%
#   ggplot(aes(x=x,y=y,z=fit,fill=fit)) +
#   geom_contour() +
#   geom_raster(alpha=0.1) +
#   facet_grid(concreteness ~ condition) +
#   scale_fill_distiller(type="div") +
#   theme_light()

#' # Session Info
#'

sessionInfo()
