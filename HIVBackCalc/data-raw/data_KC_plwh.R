############################################################
# Create a rough PLWH dataset for KC
# Based on PLWH reported for 2006-2012 in the PLoS One paper
############################################################

#############################################################
# EDIT THESE
#############################################################

# Local path of the package repo
package_repo  <- file.path(undx_repo, 'public', 'package1.0', 'HIVBackCalc')


#############################################################
# MAKE PLWH DATA
#############################################################

KCplwh <- data.frame(Year=2006:2012,
                     White=4188,
                     Black=458,
                     Hisp=572,
                     Total=5218)

#############################################################
# SAVE
#############################################################

filename1 <- file.path(package_repo, 'data-raw', 'KCplwh.csv')
filename2 <- file.path(package_repo, 'data', 'KCplwh.RData')

write.csv(KCplwh, file=filename1, row.names=FALSE)
save(KCplwh, file=filename2)


