############################################################
# This file makes a simulated KC (MSM) testing histories file 
# for HIVBackCalc function demos
#
# The code is taken from analysis_WA/combine_MSM_report.Rnw
# History: see HIVBackCalc_App/developtment/data_KC+WA_formatted.csv
############################################################

#############################################################
# EDIT THESE
#############################################################

# Local path of the undiagnosed repo
undx_repo <- '/Users/jeanette/Dropbox/School/PhD/HIV_WA'
# Local path of the package repo
package_repo  <- file.path(undx_repo, 'public', 'package1.0', 'HIVBackCalc')

#############################################################
# SETUP
#############################################################
# Package and seed
library(survival)
library(plyr)
library(Rjkb) # install using install_github
set.seed(98103)
                          
# Load data and clean data
setup_hivbackcalc(workd=undx_repo,
                  datafile='data/MSM_KingCounty_rev.csv',
                  source_these='analysis_KC/data-cleaning_JKB.R',
                  loadlib=FALSE,
                  msm=TRUE)

#############################################################
# ADDITIONAL DATA FORMATTING
#############################################################

msm <- rename(msm, c('everTested'='everHadNegTest',
                     'hiv_age_yrs'='hdx_age'))
msm$mode <- 'MSM'
msm <- transform(msm, agecat5=cut(hdx_age,
                              breaks=c(0,seq(20,70,by=5),85),
                              include.lowest=TRUE,
                              right=TRUE,
                              labels=c('<=20',
                                       '21-25',
                                       '26-30',
                                       '31-35',
                                       '36-40',
                                       '41-45',
                                       '46-50',
                                       '51-55',
                                       '56-60',
                                       '61-65',
                                       '66-70',
                                       '71-85')))

race_levels <- c('White', 'Black', 'Hisp', 'Asian', 'Native', 'Multi')
msm <- within(msm, {
                race <- as.character(racel)
                race[race=='1White'] <- 'White'
                race[race=='2Black'] <- 'Black'
                race[race=='3Hisp'] <- 'Hisp'
                race[race=='4Asian'] <- 'Asian'
                race[race=='4PI' | race == '5AmInd'] <- 'Native'
                race[race=='6Multi'] <- 'Multi'
                race <- factor(race,
                                labels=race_levels,
                                levels=race_levels)
                })
# Combine
these_cols <- c('hdx_age', 'agecat5', 'mode', 'race', 
                'everHadNegTest', 'timeDx', 'yearDx',
                'infPeriod', 'infPeriod_imputeNA')
msm <- data.frame(Population='KC', msm[,these_cols])

#############################################################
# INITIALIZE SIMULATED DATA
#############################################################
msize <- nrow(msm)

msm.sim <- data.frame(id=1:msize,
                      Population='KC-sim',
                      mode='MSM',
                      mode2='MSM')
    # 3/4/16 - added 'mode2' because we're changing the app default
    # to stratify by MSM vs non-MSM, and the default name for that
    # variable from the eHARS formatting is mode2

#############################################################
# INVESTIGATE AND ASSIGN KEY NON-TID VARIABLES
#############################################################

# Function to tabulate stuff by race
tabByRace <- function(var) {
    tabRaceVar <- table(msm$race, msm[[var]], useNA='ifany')
    tabRaceVar <- round(100*tabRaceVar/rowSums(tabRaceVar),0)
}

# Race
tabRace <- table(msm$race)
msm.sim$race <- sample(names(tabRace), size=msize,
                       prob=tabRace, replace=TRUE)

# Age
tabRaceAge <- tabByRace('agecat5')
    # From this it's a fair ROUGH approximation to say that
    # age is uniformly distributed between 21 and 45 for all
    # races
msm.sim$hdx_age <- sample(21:45, size=msize, replace=TRUE)
msm.sim <- transform(msm.sim, agecat5=cut(hdx_age,
                              breaks=c(0,seq(20,70,by=5),85),
                              include.lowest=TRUE,
                              right=TRUE,
                              labels=c('<=20',
                                       '21-25',
                                       '26-30',
                                       '31-35',
                                       '36-40',
                                       '41-45',
                                       '46-50',
                                       '51-55',
                                       '56-60',
                                       '61-65',
                                       '66-70',
                                       '71-85')))

# YearDx
tabRaceYr <- tabByRace('yearDx')
    # Same deal: it's fine to assign yearDx uniformly across
    # 2006-2012, and we'll just assign quarters randomly too
msm.sim$yearDx <- sample(2006:2012, size=msize, replace=TRUE)
msm.sim$timeDx <- msm.sim$yearDx + 
                  sample(c(0:3)/4, size=msize, replace=TRUE)

# everHadNegTest
tabRaceTest <- tabByRace('everHadNegTest')
    # Approximate this for everyone as: 
    # 10% no, 75% yes, and 15% missing
msm.sim$everHadNegTest <- sample(c(FALSE, TRUE, NA),
                                 size=msize, replace=TRUE,
                                 prob=c(0.1, 0.75, 0.15))


#############################################################
# SIMULATE infPeriod FROM RACE-SPECIFIC KM CURVES
#############################################################
# Just to be clear, I'm simulating a date of last negative
# test AND imputing the min(age-16, 17.98) rules for the 
# everHadNegTest=NOs in the variable infPeriod. This assumption
# will also be applied to the everHadNegTest=NAs in the variable
# infPeriod_imputeNA. The code will still calculate the exact
# time of infection to diagnosis using the various assumptions

# Function to use by race
sim_TID <- function(df, simdf) {

    # Subset the simdf to the race of df
    simdf <- subset(simdf, race==unique(df$race))
    
    # Subset to those who have had a negative test
    # This first step is necessary to avoid unknowingly
    # treating the NAs in everHadNegTest as TRUE
    negTest <- ifelse(df$everHadNegTest==TRUE & 
                      !is.na(df$everHadNegTest), TRUE, FALSE)
    negTestSim <- ifelse(simdf$everHadNegTest==TRUE & !  
                         is.na(simdf$everHadNegTest), TRUE, FALSE)
    time <- df$infPeriod[negTest]
    nppl <- length(time)

    # Get the survival curve
    infPeriodCurve <- summary(survfit(Surv(time, rep(1,nppl))~1))

    # Simulate from the survival curve
    # Using f=0 to take the left-sided value if the draw
    # falls between steps, arbitrarily.
    # 3/30/15: We are now also excluding simulating simtimes=0,
    # so the yright condition has been edited accordingly
    rand <- runif(n=sum(negTestSim))
    simtimes <- approx(x=infPeriodCurve$surv,
                       y=infPeriodCurve$time,
                       method='constant',
                       xout=rand,
                       f=0,
                       yleft=max(infPeriodCurve$time),
                       yright=min(infPeriodCurve$time[infPeriodCurve$time>0]))[['y']]

    # Confirm a good job
    # plot(survfit(Surv(time, rep(1,nppl))~1))
    # points(simtimes, rand, col='red')

    # Insert times back into full set of records
    simtimes.full <- rep(NA,nrow(simdf))
    simtimes.full[negTestSim] <- simtimes

    # Now insert the assumption of min(age-16, 17.98) for the
    # everHadNegTest=NOs
    noTestSim <- !negTestSim & !is.na(simdf$everHadNegTest)
    assumption <- pmin(simdf$hdx_age-16,18)
    simtimes.full[noTestSim] <- assumption[noTestSim]
    simdf$infPeriod <- simtimes.full

    # Now make infPeriod_imputeNA using that same assumption for 
    # everHadNegTest=NA
    simdf <- within(simdf, { 
                    infPeriod_imputeNA <- infPeriod
                    infPeriod_imputeNA[is.na(infPeriod)] <- 
                        assumption[is.na(infPeriod)]
                    })

    # Confirm a good job
    # table(sim=is.na(simtimes.full), everHad=is.na(simdf$everHadNegTest))
    # head(simdf[is.na(simdf$infPeriod),])

    return(simdf)
}

# Now apply the function by race
msm.sim <- ddply(msm, .(race), .fun=sim_TID, msm.sim)

# Return to original sort, just for fun
msm.sim <- msm.sim[order(msm.sim$id),]

# Rename KCsim to become compatible with code below
KCsim <- msm.sim 

###################################################################################
# Add a fake subgroup variable, for unit tests
###################################################################################
KCsim <- transform(KCsim, 
                   fakeGroup=factor(sample.int(2, nrow(KCsim), replace=TRUE),
                                    levels=c(1,2),
                                    labels=c('Group 1', 'Group 2')))

###################################################################################
# Create vectors for example with fewer repeat testers
###################################################################################

# Randomly delete prior negative tests and impute infPeriod using the standard rules.

everHadNegTest <- everHadNegTestM <- KCsim$everHadNegTest
infPeriod <- infPeriodM <- KCsim$infPeriod

table(everHadNegTest, useNA="always")
table(round(infPeriod,0), useNA="always")

for (i in 1:length(everHadNegTest)) {
  change <- rbinom(1,1,0.43)
  if (!is.na(everHadNegTest[i]) & everHadNegTest[i]==TRUE & change==1) 
    everHadNegTestM[i] <- FALSE
  
  if(!is.na(everHadNegTestM[i]) & everHadNegTestM[i]==FALSE) {
    infPeriodM[i] <- ifelse(KCsim$hdx_age[i]>33, 18, KCsim$hdx_age[i]-16)
  }
}

# check the results, should be ~= TRUE and FALSE for everHadNegTest, and
# 300+ 18 yr infPeriods.

rbind(before=table(everHadNegTest, useNA="always"), 
      after=table(everHadNegTestM, useNA="always"))
rbind(before=table(round(infPeriod,0), useNA="always"), 
      after=table(round(infPeriodM,0), useNA="always"))

# add new vectors to dataframe

KCsim <- cbind(KCsim, everHadNegTestM, infPeriodM)

#############################################################
# SAVE
#############################################################

filename1 <- file.path(package_repo, 'data-raw', 'KCsim.csv')
filename2 <- file.path(package_repo, 'data', 'KCsim.RData')

write.csv(KCsim, file=filename1, row.names=FALSE)
save(KCsim, file=filename2)

