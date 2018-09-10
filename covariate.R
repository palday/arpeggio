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
         # you can insert more meaningful labels for the conditions here
         condition = factor(condition, levels=c(1,2,3), labels=c("A","B","C")))

#' Load and convert channel coordinates from spherical to cartesian.
#' see https://sccn.ucsd.edu/pipermail/eeglablist/2006/001655.html for formula
channel_locations <- read.table("64Ch_actiCAP_Equidistant10_polar.txt",header=TRUE,sep="\t")

#' In Cartesian coordinates, (0,0,0) is the center of the skull and (0,0,1) is
#' the apex, i.e. electrode Cz in the 10-20 system. Negative x are left lateral,
#' positive x are right lateral. Positive y is anterior, negative y is
#' posterior.

polar2cart <- function(phi,theta, r=1){
  list(x = r * sin(phi) * cos(theta),
       y = r * sin(phi) * sin(theta),
       z = r * cos(phi))
}

channel_locations <- channel_locations %>%
  do(as.data.frame(polar2cart(.$phi,.$theta))) %>%
  bind_cols(channel_locations) %>%
  mutate(chan = as.character(chan))

#' Now add in topographical information to our original data
dat <- dat %>%
  # strip out the C prefix before the channels to match our coordinates
  mutate(chan=substring(chan,2)) %>%
  left_join(channel_locations)

#' Item number is missing here. Items aren't the same as trials -- the same
#' items can be presented in different orders, and it's not clear if trlnumber
#' here refers to the sequence number or the particular lexical item. Items
#' should definitely be included in the random effects the same way subjects
#' are. Although given the number of controls implemented in the covariates, you
#' should be able to get away with (1|item) instead of (1+condition|item).

m <- lmer(N4 ~ BS + condition * x * y *
                  (cloze  + wordfrq + phon_nd + pt_semdist + sentt_semdist + concreteness + rhyme_eval + plaus_eval) +
              (1 + condition | subject_id),
          data=dat)

summary(m)

Anova(m)

#' For now, I'm leaving out ways of plotting the statistical model, but that's
#' something we can add later. Including potentially ways to plot the
#' coefficients as ERP topographies.
