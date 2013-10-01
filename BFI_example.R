#########################################################################################
#
# Example for Baseflow seperation using the Lyne and Hollick filter
#
# tony.ladson@gmail.com
# 1 Oct 2013
#
########################################################################################
#
# An example using the BFI function at 
# https://github.com/TonyLadson/BaseflowSeparation_LyneHollick

library(RCurl)
library(foreign)

# source the BFI function
source('https://raw.github.com/TonyLadson/BaseflowSeparation_LyneHollick/master/BFI.R')

args(BFI)
# function (Q, alpha = 0.925, passes = 3, ReturnQbase = FALSE, 
#           n.reflect = 30) 


# Get the sample data set
# 67 daily flow values for the Bass River at Loch
# This is the same data set provided in the paper
# Ladson, A. R., R. Brown, B. Neal and R. Nathan (2013)
# A standard approach to baseflow separation using the Lyne and Hollick filter. Australian Journal of Water Resources 17(1): 173-180.
# http://dx.doi.org/10.7158/W12-028.2013.17.1

url <- "https://raw.github.com/TonyLadson/data/master/data/BassRiver@Loch.csv"

# I needed to tell R where to find the certificates to download data from github
ca.path <- "C:/Users/TonyLadson/Documents/R/win-library/3.0/RCurl/CurlSSL/cacert.pem"

# Q is a vector of flow values for the Bass River at Loch
Q <- as.vector(unlist(read.csv(textConnection(getURL(url, cainfo=ca.path)), header=FALSE)))

plot(Q, type='l')

BFI(Q, alpha=0.98)

# $BFI
# [1] 0.1965759
# 
# $alpha
# [1] 0.98
# 
# $FractionUsed
# [1] 1

lines(BFI(Q, alpha=0.98, ReturnQbase=TRUE)$Qbase, lty=2)
lines(BFI(Q, alpha=0.925, ReturnQbase=TRUE)$Qbase, lty=2, col=4)


########## Example with Missing Values ##############


Q <- c(NA,5,7,108,117,57,36,26,95,1169,308,
       144,89,62,48,40,35,73,82,342,393,310,
       275,260,245,256,141,119,934,382,158,96,
       122,103,83,67,148,NA,366,161,119,82,330,294,
       261,266,153,247,703,498,286,163,124,85,94,
       81,62,47,37,30,26,24,24,22,21,20,19,18,18,17,16,NA,20,19,18,18,NA)

BFI(Q, alpha=0.98)

$BFI
# [1] 0.1411197
# 
# $alpha
# [1] 0.98
# 
# $FractionUsed
# [1] 0.9452055


# Explanation of FractionUsed

#   Number of daily values
sum(!is.na(Q))
#    73

#   with n.reflect set to 30 (the default), the final 4 values between the NA can not be used
#   therefore the faction used is
(73-4)/73

# The BFI is the weighted average of the BFI of the flow segments with non-missing values

(length(Q[2:37])*BFI(Q[2:37], alpha=0.98)$BFI + length(Q[39:71])*BFI(Q[39:71], alpha=0.98)$BFI)/(length(Q[2:37]) + length(Q[39:71]))
# [1] 0.1411197

########## Example with more passes ##############


url <- "https://raw.github.com/TonyLadson/data/master/data/BassRiver@Loch.csv"
ca.path <- "C:/Users/TonyLadson/Documents/R/win-library/3.0/RCurl/CurlSSL/cacert.pem"
Q <- as.vector(unlist(read.csv(textConnection(getURL(url, cainfo=ca.path)), header=FALSE)))

plot(Q, type='l')
lines(BFI(Q, ReturnQbase=TRUE, passes=3)$Qbase, lty=2)
lines(BFI(Q, ReturnQbase=TRUE, passes=9)$Qbase, lty=2, col=4)