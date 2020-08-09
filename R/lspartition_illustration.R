################################################################################
# lspartition: illustration file
# Authors: M. D. Cattaneo, M. Farrell and Y. Feng
# Last update: 07-AUG-2019
################################################################################
rm(list=ls(all=TRUE))
library(TeachingDemos); library(lspartition); library(ggplot2)
set.seed(1234)

#######################################################
## SET =TRUE TO GENERATE OUTPUTS
do.output = TRUE
#######################################################

## Install NPPACKAGE Package
if(do.output) txtStart("output/lspartition_1.txt", results = FALSE)
install.packages("lspartition", dependencies = TRUE)
library(lspartition)
if(do.output) txtStop()

## read data
if(do.output) txtStart("output/lspartition_2.txt")
data <- read.csv("bikesharing.csv", header = TRUE)
summary(data)
if(do.output) txtStop()

## outcome, covariate
if(do.output) txtStart("output/lspartition_3.txt")
y <- data$count
x <- data$atemp
g <- data$workingday
if(do.output) txtStop()

####kappa selection
if(do.output) txtStart("output/lspartition_4.txt")
summary(lspkselect(y, x, kselect="imse-rot", subset=(g==1)))
if(do.output) txtStop()

if(do.output) txtStart("output/lspartition_5.txt")
summary(lspkselect(y, x, kselect="imse-dpi", subset=(g==1)))
if(do.output) txtStop()

if(do.output) txtStart("output/lspartition_6.txt", result = FALSE)
summary(lspkselect(y, x, kselect="imse-dpi", ktype="qua", subset=(g==1)))
if(do.output) txtStop()

## pointwise inference
if(do.output) txtStart("output/lspartition_7.txt")
est_workday_bc1 <- lsprobust(y, x, neval=20, bc= "bc1", nknot=8, subset = (g==1))
est_workday_bc3 <- lsprobust(y, x, neval=20, bc= "bc3", nknot=8, bnknot=10, subset = (g==1))
summary(est_workday_bc1)
if(do.output) txtStop()

if(do.output) txtStart("output/lspartition_8.txt")
lsprobust.plot(est_workday_bc1, xlabel="Temperature", ylabel="Number of Rentals", legendGroups = "Working Days") + theme(text=element_text(size=17), legend.position=c(.15,.9))
ggsave("output/pointwise1.pdf", width=6.8, height=5.5)
lsprobust.plot(est_workday_bc3, xlabel="Temperature", ylabel="Number of Rentals") + theme(text=element_text(size=17), legend.position="none")
ggsave("output/pointwise2.pdf", width=6.8, height=5.5)
if(do.output) txtStop()


## uniform inference: numerator matrix
if(do.output) txtStart("output/lspartition_9.txt")
est_workday_bc1 <- lsprobust(y, x, bc = "bc1", nknot = 4, uni.method = "pl", uni.ngrid = 100, uni.out = T, subset = (g==1))
round(est_workday_bc1$uni.output$t.num.pl[1:5,],3)
if(do.output) txtStop()

## uniform inference: plug-in method
if(do.output) txtStart("output/lspartition_10.txt")
est_workday_bc1 <- lsprobust(y, x, neval=20, bc= "bc1", uni.method="pl", nknot=8, subset = (g==1), band = T)
est_workday_bc1$sup.cval
if(do.output) txtStop()

if(do.output) txtStart("output/lspartition_11.txt")
lsprobust.plot(est_workday_bc1, CS="all", xlabel="Temperature", ylabel="Number of Rentals", legendGroups = "Working Days") + theme(text=element_text(size=17), legend.position=c(.15,.9))
ggsave("output/uniform1.pdf", width=6.8, height=5.5)
if(do.output) txtStop()

## bootstrap
if(do.output) txtStart("output/lspartition_12.txt")
est_workday_bc3 <- lsprobust(y, x, neval=20, bc= "bc3", nknot=8, bnknot=10, uni.method="wb", subset = (g==1), band = T)
est_workday_bc3$sup.cval
lsprobust.plot(est_workday_bc3, CS="all", xlabel="Temperature", ylabel="Number of Rentals", legendGroups = "Working Days") + theme(text=element_text(size=17), legend.position=c(.15,.9))
ggsave("output/uniform2.pdf", width=6.8, height=5.5)
if(do.output) txtStop()


####Two groups: pointwise
if(do.output) txtStart("output/lspartition_13.txt")
est_workday  <- lsprobust(y, x, neval=20, bc= "bc3", nknot=8, subset = (g==1))
est_nworkday <- lsprobust(y, x, neval=20, bc= "bc3", nknot=8, subset = (g==0))
lsprobust.plot(est_workday, est_nworkday, legendGroups=c("Working Days", "Nonworking Days"), xlabel="Temperature", ylabel="Number of Rentals", lty=c(1,2)) + theme(text=element_text(size=17), legend.position=c(.2,0.85))
ggsave("output/diff1.pdf", width=6.8, height=5.5)
if(do.output) txtStop()

## Two groups: diff
if(do.output) txtStart("output/lspartition_14.txt")
diff <- lsplincom(y, x, data$workingday, R=c(-1,1), band=T, cb.method="pl")
summary(diff)
if(do.output) txtStop()

if(do.output) txtStart("output/lspartition_15.txt")
lsprobust.plot(diff, CS="all", xlabel="Temperature", ylabel="Number of Rentals", legendGroups="Difference between Working and Other Days") + theme(text=element_text(size=17), legend.position=c(.36,.2))
ggsave("output/diff2.pdf", width=6.8, height=5.5)
if(do.output) txtStop()

## Two groups: diff, smoother fit
if(do.output) txtStart("output/lspartition_16.txt")
diff <- lsplincom(y, x, data$workingday, R=c(-1,1), band=T, cb.method="pl", m=3)
lsprobust.plot(diff, CS="all", xlabel="Temperature", ylabel="Number of Rentals") + theme(text=element_text(size=17), legend.position="none")
ggsave("output/diff3.pdf", width=6.8, height=5.5)
if(do.output) txtStop()

###################################################################
# The following section shows how to manually depict the estimated   
# curve, confidence region, and true function in the same plot. 
# It is NOT reported in the paper. 
###################################################################

# Linear spline fit using lsprobust()
est_workday_bc1 <- lsprobust(y, x, neval=20, bc= "bc1", uni.method="pl", nknot=8, subset = (g==1), band = T)
xeval   <- est_workday_bc1$Estimate[,"X1"]           # evaluation points
yhat    <- est_workday_bc1$Estimate[,"tau.cl"]       # fitted value, uncorrected
yhat.bc <- est_workday_bc1$Estimate[,"tau.bc"]       # fitted value, bias corrected
se      <- est_workday_bc1$Estimate[,"se.rb"]        # standard errors, bias corrected
cval    <- est_workday_bc1$sup.cval                  # critical value for confidence band

# A global polynomial fit, treated as the true function
polyfit <- lm(y~x+I(x^2)+I(x^3)+I(x^4), subset=(g==1))$coefficients
y.true  <- polyfit[1]+xeval*polyfit[2]+xeval^2*polyfit[3]+xeval^3*polyfit[4]+xeval^4*polyfit[5]

##############################
# Plotting using R base graphs
# True function
plot(xeval, y.true, type="l", col="purple")
# fitted values
lines(xeval, yhat, lty=3)
# confidence intervals
segments(xeval, yhat.bc-1.96*se, xeval, yhat.bc+1.96*se, col="darkgrey")
# confidence band
lines(xeval, yhat.bc-cval*se, lty=2, col="navy")
lines(xeval, yhat.bc+cval*se, lty=2, col="navy")

#############################
# Alternatively, adding the true function after lsprobust.plot()
fig <- lsprobust.plot(est_workday_bc1, CS="all", xlabel="Temperature", ylabel="Number of Rentals")+theme(text=element_text(size=17), legend.position="none")
fig <- fig +geom_line(data=data.frame(x=xeval, y=y.true), aes(x=xeval, y=y.true), col="purple", lty=2)
