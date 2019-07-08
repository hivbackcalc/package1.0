
#' Aggregates diagnosis counts into the specified time interval
#'  
#' HIV diagnosis counts from timeDx are aggregated to match 
#' the specified intLength. Called by tabulateDiagnoses()
#'  
#' @param timeDx Time of diagnosis from the testhist data. 
#'        Form should be either YYYY+x/12 where x is the month number, YYYY.25,
#'        YYYY.5, or YYYY. 
#' @param intLength Desired interval length for diagnoses: 1/12, 
#'        0.25, 0.5 or 1 (1=1 year)
#'  
#' @return timeDx variable with intLength 
#' @seealso \code{\link{tabulateDiagnoses}} for example

aggregateDiagnoses <- function(timeDx, intLength) {

    # Don't allow irregular time steps
    if (!intLength%in%c(1/12, 0.25,0.5,1)) stop('Error in tabulateDiagnoses:
                                          intLength must be 1/12, 0.25, 0.5 or 1')


    
    # Check that intLength is >= time steps in timeDx
    times <- sort(unique(timeDx))
    diffs <- diff(times)
    tstep <- median(diffs)
    check <- round(intLength-tstep, 4)
    if (check<0) warning('Error in tabulateDiagnoses: intLength is smaller
                              than median interval between diagnosis counts
                              in variable timeDx')
    
    # Note the year
    floored <- floor(timeDx)

    # Either leave timeDx alone, or aggregate to the desired intLength
    if (abs(check)<.0001) {
        timeDxNew <- timeDx
    } else if (intLength %in% c(0.25, 0.5, 1)) {
        warning('\n Because intLength=', intLength, ', aggregating diagnosis times to ', intLength, '-years')
        decimals <- timeDx-floored
        timeDxNew <- floored + floor(decimals/intLength)*intLength 
    } else stop('Error in aggregateDiagnoses(): timeDx and/or intLength incorrectly specified')

    return(timeDxNew)
}


#' Format eHARS person view data for use with the HIVBackCalc package
#'  
#' Creates a data frame of class 'testinghistories' with formatted and
#' calculated variables. Additionally returns a record of formatting
#' assumptions made for inconsistent entries
#'  
#' @param rawdata Unformatted data frame or file path to CSV file
#' @param assumptionNo Choice of assumption for those reporting never 
#'          having had a negative test. Default is 'age16', which imputes
#'          time from infection to diagnosis as min(age-16, 18 yrs). 'age16mid'
#'          instead uses the midpoint between age and 16, i.e. 
#'          min((age-16)/2, 18 yrs). 
#' @return Data frame of class "testinghistories" 
format_eHARS <- function(rawdata, assumptionNo='age16') {
  
    require(plyr)

    # Read in data
    if (!is.data.frame(rawdata)) {
        rawdata  <- read.csv(rawdata,
                             na.strings='',
                             stringsAsFactors=FALSE)
    }

    # Helper function to flag records for cleaning
    recordFlag  <- function(df, logical, message) {
        # Set NA's in logical to FALSE
        logical[is.na(logical)] <- FALSE
        # Create flag variable if it doesn't exist
        if (!'flag'%in%colnames(df)) df$flag=''
        # Record flag and tidy up
        df$flag[logical] <- paste0(df$flag[logical], '; ', message) 
        df$flag <- gsub('^;', '', df$flag)
        # Add to assumptions table
        assumptionsN[assumptionsN$Issue==message,'N'] <<- sum(logical) 
        return(df)
    }

    # Construct empty assumptions table
    assumptionsN <- data.frame(Issue=c('Missing month',
                                    'Missing day',
                                    'Illogical last negative',
                                    'everHadNegTest inconsistent with infPeriod',
                                    'infPeriod capped at aidsUB',
                                    'Age <=16 and no infPeriod'),
                               Assumption=c('Month (diagnosis or last neg test) assumed July for computing infPeriod; diagnosis quarter randomly imputed',
                                            'Day (diagnosis or last neg test) assumed 15th for computing infPeriod',
                                            'Last negative date overwritten as missing because recorded as at or after diagnosis',
                                            'everHadNegTest flag altered to match presence/absence of last neg test date',
                                            'infPeriod capped at 18 years',
                                            'Removed from dataset because <=16 yrs and no recorded infPeriod'),
                               N=NA)

    #############################################################
    # OVERVIEW: N, expected variables and labels for factors
    #############################################################
    dataf <- rawdata
    colnames(dataf) <- tolower(colnames(dataf))
    variables_expected <- c('rsh_state_cd',
                          'aids_dx_dt',
                          'cd4_first_hiv_dt',
                          'cd4_first_hiv_type',
                          'cd4_first_hiv_value',
                          'hiv_aids_age_yrs',
                          'hiv_dx_dt',
                          'race',
                          'screen_last_neg_dt',
                          'trans_categ',
                          'vl_first_det_dt',
                          'vl_first_det_value',
                          'birth_sex',
                          'stateno',
                          'tth_last_neg_dt',
                          'tth_first_pos_dt',
                          'tth_ever_neg')
    not_in_dataf <- variables_expected[!variables_expected%in%colnames(dataf)]
    if (length(not_in_dataf)!=0) stop('Some eHARS variables are missing: \n', 
                                      paste(not_in_dataf,
                                            collapse='\n'))
    
    # Factors and labels
    races <- c('Hispanic',
                'American Indian/Alaska Native',
                'Asian',
                'Black',
                'Native Hawaiian/Pacific Islander',
                'White',
                'Legacy Asian/Pacific Islander',
                'Multi-race',
                'Unknown')
    modes <- c('Adult MSM', 
                'Adult IDU',
                'Adult MSM & IDU',
                'Adult received clotting factor',
                'Adult heterosexual contact',
                'Adult received transfusion/transplant',
                'Perinatal exposure, HIV diagnosed at age 13 years or older',
                'Adult with other confirmed risk',
                'Adult with no identified risk (NIR)',
                'Adult with no reported risk (NRR)',
                'Child received clotting factor',
                'Perinatal exposure',
                'Child received transfusion/transplant',
                'Child with other confirmed risk',
                'Child with no identified risk (NIR)',
                'Child with no reported risk (NRR)',
                'Risk factors selected with no age at diagnosis')
    dataf <- transform(dataf,
                     new_race=as.character(factor(race, levels=1:9,
                                     labels=races)),
                     new_mode=as.character(factor(trans_categ, 
                                                  levels=c(1:13,18:20,99), 
                                                  labels=modes)),
                       stringsAsFactors=FALSE)
    # Some renaming
    dataf <- rename(dataf,c('hiv_aids_age_yrs'='hdx_age',
                              'birth_sex'='sex'))
    

    #############################################################
    # SUMMARIZE: Summarize or tabulate variables and store in
    #            a table
    #############################################################
    varsummaries <- vector('list', length=ncol(dataf))
        names(varsummaries) <- colnames(dataf)
        for (x in 1:ncol(dataf)) {
        varname <- colnames(dataf)[x]
        var = dataf[,x]
        if (length(unique(var))>25) {
          if (is.numeric(var)) {
            sum <- summary(var)            
            if ("NA's"%in%names(sum)) {
                nmiss <- sum["NA's"]
            } else nmiss <- 0
            sum <- sum[c('Min.', 'Mean', 'Max.')]
          } else {
            nmiss <- table(var)[is.na(names(table(var)))]
            sum <- ''
          }
        } else {
          sum <- table(var, useNA='always')
          nmiss <- sum[is.na(names(sum))]
        }
        nmiss <- as.numeric(nmiss)
        result <- c(sum, round(100*nmiss/nrow(dataf),2))
        names(result)[length(result)] <- 'Percent Miss'
        result <- result[!is.na(names(result))]
        varsummaries[[x]] <- data.frame(Variable=c(varname, 
                                                   rep('', length(result)-1)),
                                        Values=names(result),
                                        N=result,
                                        row.names=NULL)
    }
    varsummaries <- data.frame(do.call('rbind', varsummaries),
                               row.names=NULL)


    #############################################################
    # COLLAPSE RACE AND MODE OF DIAGNOSIS
    #############################################################

    collapsed_race <- c('White', 'Black', 'Hispanic', 'Asian', 
                        'Native', 'Multi/Other')
    collapsed_mode <- c('MSM', 'Hetero', 'Blood/Needle/Other')
    dataf <- within(dataf, {
        race6 <- as.character(new_race)
        race6[race6 %in% c("American Indian/Alaska Native", 
                           "Native Hawaiian/Pacific Islander")] <- 'Native'
        race6[race6 %in% "Legacy Asian/Pacific Islander"] <- 'Asian'
        race6[race6 %in% c("Multi-race","Unknown")] <- 'Multi/Other'
        mode3 <- as.character(new_mode)
        mode3[mode3 %in% c('Adult MSM','Adult MSM & IDU')] <- 'MSM'
        mode3[mode3 %in% c('Adult heterosexual contact',
                           'Adult with no identified risk (NIR)')] <- 'Hetero'
        mode3[!mode3 %in% c('MSM', 'Hetero')] <- 'Blood/Needle/Other'
#        race6 <- factor(race6,
#                       labels=collapsed_race,
#                       levels=collapsed_race)
#        mode3 <- factor(mode3,
#                       levels=collapsed_mode,
#                       labels=collapsed_mode)
#        mode2 <- factor(ifelse(mode3 %in% 'MSM', 'MSM', 'non-MSM'))
         mode2 <- ifelse(mode3 %in% 'MSM', 'MSM', 'non-MSM')
        # FOR NOW: make the main mode=mode2
        mode <- mode2
    })

    #############################################################
    # DEFINE AGE GROUPS
    #############################################################
    dataf <- transform(dataf,
                     agecat5=cut(hdx_age,
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

    #############################################################
    # FORMAT DATES AND CREATE INFPERIOD
    #############################################################

    # Helper function to work with dates
    # For each date, need a fake date if month and/or day are missing 
    # plus an imputed quarter if month is missing

    get_dates <- function(timevar) {
        year <- suppressWarnings(as.numeric(substr(timevar,1,4)))
        month <- substr(timevar,5,6)
        day <- substr(timevar,7,8)
        missing_month <- !is.na(year) & month=='..'
        missing_day <- !is.na(year) & day=='..' & !month=='..'
        # Create a year-quarter variable, imputing a quarter if necessary
        set.seed(98103)
        yrqtr <- year + 
            suppressWarnings(ifelse(missing_month, sample(c(0,0.25,0.5,0.75)), 
                                    floor(as.numeric(month)/4)*0.25))
        # Create an  _imputed date for calculating inter-test intervals
        # 15th of the month if only month is known; July 1 if only year known
        day <- ifelse(missing_month, '01', ifelse(missing_day, '15', day))
        month <- ifelse(missing_month, '07', month)
        dateChar <- apply(cbind(year,month,day),1,paste,collapse='-')
        dateChar[dateChar=='NA-NA-NA'] <- ''
        dateImp <- as.Date(dateChar,"%Y-%m-%d")
        return(list(dateImp=dateImp, 
                    year=year,
                    yrqtr=yrqtr,
                    missMonth=missing_month, 
                    missDay=missing_day))
    }

    # Diagnosis date
    dxDate <- get_dates(dataf$hiv_dx_dt)
    # Last negative test date
    negDate <- get_dates(dataf$tth_last_neg_dt)
    dataf <- transform(dataf,
                       yearDx=dxDate$year,
                       timeDx=dxDate$yrqtr,
                       infPeriod=as.numeric(dxDate$dateImp-
                                            negDate$dateImp)/365,
                       stringsAsFactors=FALSE)
    # Record assumptions
    dataf <- recordFlag(dataf, dxDate$missMonth | negDate$missMonth,
                        'Missing month')
    dataf <- recordFlag(dataf, dxDate$missDay | negDate$missDay,
                        'Missing day')

    # Illogical last negative
    dataf <- recordFlag(dataf, dataf$infPeriod<=0,
                        'Illogical last negative')
    dataf <- within(dataf, {
                    infPeriod[infPeriod<=0] <- NA
                       })


    #############################################################
    # CREATE everHadNegTest after creating infPeriod
    #############################################################
    # Define everHadNegTest based on tth_ever_neg
    dataf <- transform(dataf, 
                     everHadNegTest=ifelse(tth_ever_neg=='Y', TRUE, 
                                           ifelse(tth_ever_neg=='N', FALSE, NA)))
    #with(dataf,table(everHadNegTest, tth_ever_neg, useNA='always'))

    # Look at actual infPeriod values by everHadNegTest
    #ddply(dataf, .(everHadNegTest), function(x) c(summary(x$infPeriod)))

    ## ---- fix_everHadNegTest_toTRUE ----
    toTRUE1 <- !dataf$everHadNegTest & !is.na(dataf$infPeriod)
    toTRUE2 <- is.na(dataf$everHadNegTest) & !is.na(dataf$infPeriod)
    dataf$everHadNegTest[toTRUE1] <- TRUE
    dataf$everHadNegTest[toTRUE2] <- TRUE

    ## ---- fix_everHadNegTest_toFALSE ----
    toFALSE <- dataf$everHadNegTest & is.na(dataf$infPeriod)
    dataf$everHadNegTest[toFALSE] <- FALSE

    # Record assumptions and flag
    dataf <- recordFlag(dataf, toTRUE1 | toTRUE2 | toFALSE,
                     'everHadNegTest inconsistent with infPeriod')
              
    ## ---- check_everHadNegTest ----
    #checkEver <- with(dataf,table(everHadNegTest, 
    #                             TID_NA=is.na(infPeriod), useNA='always')))

    #############################################################
    # EDIT infPeriod
    #############################################################

    # Cap at AIDS upper bound of ~18 years
    aidsUB <- qweibull(.95,shape=2.516,scale=1/0.086) #17.98418
    infTemp <- dataf$infPeriod
    
    fixNo <- function(dataf, assumptionNo) {
        switch(assumptionNo,
               age16 = {
                    dataf <- transform(dataf,
                                       infPeriod=ifelse(everHadNegTest, 
                                                        pmin(infPeriod, aidsUB), 
                                                        ifelse(!everHadNegTest, 
                                                               pmin(hdx_age-16, 
                                                                    aidsUB), 
                                                               NA)))
               },
               age16mid = {
                    dataf <- transform(dataf,
                                       infPeriod=ifelse(everHadNegTest, 
                                                        pmin(infPeriod, aidsUB), 
                                                        ifelse(!everHadNegTest, 
                                                               pmin((hdx_age-16)/2, 
                                                                    aidsUB), 
                                                               NA)))
               }
               )
        return(dataf)
    }
    dataf <- fixNo(dataf, assumptionNo)

    dataf <- recordFlag(dataf, infTemp!=dataf$infPeriod,
                        'infPeriod capped at aidsUB')

    # Remove cases who are too young for the impute-infPeriod assumption
    dataf <- recordFlag(dataf, 
                        dataf$hdx_age<=16 & (!dataf$everHadNegTest | 
                                        is.na(dataf$everHadNegTest)),
                        message='Age <=16 and no infPeriod')
    dataf <- subset(dataf, !(hdx_age<=16 & (!everHadNegTest | 
                                            is.na(everHadNegTest))))

    #############################################################
    # CREATE infPeriod_imputeNA
    #############################################################
    dataf <- within(dataf,{ 
        infPeriod_imputeNA=ifelse(is.na(everHadNegTest),
                                  pmin(hdx_age-16, aidsUB),
                                  infPeriod)
    })


    class(dataf) <- append(class(dataf), 'testinghistories')
    return(list(data=dataf,
                assumptions=assumptionsN,
                rawVarSum=varsummaries))
}


#' Plots diagnosis counts over time
#'
#' Plot of diagnoses vs time, with option to panel by subgroups
#'
#' @param testhist Data frame of class 'testinghistories' with variable
#'        timeDx 
#' @param showlegend Logical to turn on the legend
#' @param legendpos Character specifying legend position: bottom, top, right
#' @examples
#' plotDiagnoses(KCsim)
#' plotDiagnoses(KCsim, showlegend=TRUE) + aes(color=race)
#' plotDiagnoses(KCsim, showlegend=TRUE) + aes(color=race) + facet_grid(.~everHadNegTest)

plotDiagnoses  <- function(testhist, showlegend=FALSE, legendpos='bottom') {

    # Could use geom_count() instead...
    testhist$count <- 1
    p <- ggplot(testhist, aes(x=timeDx, y=count, color=factor(count))) +
        stat_summary(fun.y=sum, geom='point') +
        stat_summary(fun.y=sum, geom='line') +
        theme_bw() + 
        scale_color_discrete(name="") +
        theme(axis.text.x=element_text(angle=90)) + 
        scale_x_continuous(breaks=seq(min(testhist$yearDx),max(testhist$yearDx),by=1))+
        xlab("Time") + ylab("Diagnoses") 
    if (showlegend) {
        p <- p + theme(legend.position=legendpos)
    } else p <- p + guides(color=FALSE)

    return(p)
  }


#' Tabulate responses to 'Have you ever had a negative test?'
#'
#' Tabulates the everHadNegTest variable by any number of stratification/
#' subgroup variables
#'
#' @param testhist Data frame of class 'testinghistories' with variable
#'        everHadNegTest 
#' @param variables NAMED character vector of stratification variable names. 
#'        Names will be used as table labels.
#'        Variables must exist in the testhist data frame.
#'        See examples for using a list structure to do multi-variale cross-tabs.
#' @param supercolumn Set to TRUE to include a pretty column indicating
#'        stratification variables
#' @param fullsample_row Set to TRUE to have the 1st row be the results
#'        for the entire sample
#' @examples
#' # One variable
#' tabTestHist(KCsim, c(Year='yearDx'))
#' # Marginal, two variables
#' tabTestHist(KCsim, c(Year='yearDx', Race='race'), supercolumn=TRUE)
#' # Cross-tab, two variables
#' tabTestHist(KCsim, list(c(Year='yearDx', Race='race')))

tabTestHist <- function(testhist, variables, 
                        supercolumn=FALSE, fullsample_row=FALSE) {

    # variables is either a named vector or a list of length 1 of a 
    # named vector
    if (!is.list(variables)) {
        checkvec <- variables
    } else if (is.list(variables)) {
        checkvec <-  variables[[1]]
    } else stop('In tabTestHist, variables must be list or vector')
    if (is.null(names(checkvec))) stop('In tabTestHist, Vector elements must be named')
    
    # supercolumn must be TRUE if variables is a vector with >1 element
    if (is.vector(variables) & length(variables)>2) supercolumn <- TRUE
    vars <- list(NULL)
    for (v in 1:length(variables)) {

    # Tabulate everHadNegTest for this subgroup
    tab <- ddply(testhist, variables[[v]], function(x, TN=nrow(testhist)) {
        n <- nrow(x)
        c(N=n,
        `Column Percent`=round(100*n/TN,0),
        `Percent Yes`=round(100*sum(x$everHadNegTest, na.rm=TRUE)/n,0),
        `Percent No`=round(100*sum(!x$everHadNegTest, na.rm=TRUE)/n,0),
        `Percent Missing`=round(100*sum(is.na(x$everHadNegTest))/n,0))
    })
    # Add pretty column with subgroup name
    if (supercolumn) {
    colnames(tab)[1] <- 'Value'
    tab <- data.frame(Variable=rep('',nrow(tab)),
                    tab,
                    stringsAsFactors=FALSE,
                    check.names=FALSE)
    tab$Variable[1] <- names(variables)[v]
    }

    vars[[v]] <- tab
  }
  
  # Compile results into one table
  vars <- do.call(rbind, vars)
  
  # Now add a row for the full sample
  if (fullsample_row) {
    fulleverHadNegTest <- round(100*table(testhist$everHadNegTest, 
                                          useNA='ifany')/nrow(testhist),0)
    fullrow <- data.frame(vars[1,])
    thiscol <- which(colnames(fullrow)=='N')
    fullrow[,c((ncol(fullrow)-4):ncol(fullrow))] <-
      c(nrow(testhist), 
        100, 
        fulleverHadNegTest['TRUE'], 
        fulleverHadNegTest['FALSE'],
        100-fulleverHadNegTest['TRUE']-fulleverHadNegTest['FALSE'])
    fullrow[,1:(thiscol-1)] <- rep('All', thiscol-1)
    colnames(fullrow) <- colnames(vars)
    vars <- rbind(fullrow, vars)
  }
  
  return(vars)
}


#' Plot responses to 'Have you ever had a negative test?' over time
#'
#' Plots the everHadNegTest results over time
#'
#' @param testhist Data frame of class 'testinghistories' with variables
#'        everHadNegTest and yearDx
#' @param group Optional grouping variable to display using different line types
#' @param panel Optional variable by which to panel the plot
#' @examples
#' # Simple
#' plotTestHist(KCsim)
#' # Add a simple group. Refer to it as Group when calling the plot.
#' newdat <- transform(KCsim, fakeGroup=factor(sample.int(2, nrow(KCsim), replace=TRUE)))
#' plotTestHist(newdat, group='fakeGroup') + aes(lty=Group)
#' # That's messy - let's do panels instead. Refer to it as Panel when calling the plot
#' plotTestHist(newdat, panel='fakeGroup') + facet_grid(.~Panel)

plotTestHist <- function(testhist, group=NULL, panel=NULL) {

    id <- c(Year='yearDx', Group=group, Panel=panel)
    id <- id[!is.null(id)]

    # Tabulate testing history breakdown
    tabTime <- tabTestHist(testhist, list(id))

    # Drop the column percent and N variables, and melt
    tabTime <- subset(tabTime, select=-c(`Column Percent`, N))
    colnames(tabTime) <- gsub('Percent ','',colnames(tabTime))
    tabTime <- melt(tabTime, id.vars=names(id))

    p <- ggplot(tabTime,aes(x=Year,y=value,color=variable))  +   
        geom_point() + geom_line() +
        theme_bw()+
        theme(legend.position='bottom',axis.text.x=element_text(angle=90)) + 
        scale_color_hue(name="Ever had negative test?") + 
        scale_x_continuous(breaks=seq(min(tabTime$Year),max(tabTime$Year),by=2))+
        xlab("Time") + ylab("Percent") 

    return(p)
}

