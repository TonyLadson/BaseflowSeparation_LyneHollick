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

