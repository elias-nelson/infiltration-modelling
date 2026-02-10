################################################################################
#
#  InfiltraStat: R functions for infiltration modelling, geostatistics,
#  and spatial soil analysis developed for research on soil water dynamics
#  and conservation agriculture systems.
#
#  Author: Elias Nelson Kaunda
#  Inspired in part by methodologies from the CEPHaStat / GeoNutrition
#  projects (Murray Lark et al.), released under the GNU GPLv3.
#
#  This collection provides tools for:
#    - Statistical analysis of soil hydraulic and environmental data
#    - Infiltration model fitting and evaluation
#    - Variogram modelling and spatial prediction
#    - Mapping of infiltration and kriging uncertainty
#
#  License: GNU General Public License v3.0
#
################################################################################

vers <- "v0.1 Initial development:"

header <- "
InfiltraStat is a set of R functions developed for analysing soil
infiltration processes, hydraulic behaviour, and spatial variability
in agricultural landscapes. The code is distributed under the terms
of the GNU General Public License v3.0, without liabilty and warranty.
"

fnclist <- "
Current function groups:
  - Descriptive statistics and diagnostics
  - Infiltration model fitting (Horton, Philip, Green-Ampt)
  - Variogram estimation and kriging
  - Spatial mapping and uncertainty analysis
"


cat(paste("InfiltraStat version ",vers))
cat(header)
cat(fnclist)

#-------------------------------------------------------------------------------
# Compute sample skewness of a numeric vector
#
# Uses bias-corrected moment coefficient of skewness.
# Normal distribution -> skewness = 0
#-------------------------------------------------------------------------------
skew <- function(x, na.rm = TRUE) {
  
  if (na.rm) x <- stats::na.omit(x)
  
  n <- length(x)
  
  if (n < 3)
    stop("skew() requires at least 3 observations.")
  
  mu <- mean(x)
  sd_x <- stats::sd(x)
  
  if (sd_x == 0)
    stop("skew undefined for zero-variance data.")
  
  sk <- (n / ((n - 1) * (n - 2))) *
    sum(((x - mu) / sd_x)^3)
  
  return(sk)
}

#-------------------------------------------------------------------------------
# Compute kurtosis of a numeric vector
#
# Returns excess kurtosis by default.
#
# Normal distribution:
#   - Pearson kurtosis = 3
#   - Excess kurtosis  = 0
#-------------------------------------------------------------------------------
kurt <- function(x, na.rm = TRUE, excess = TRUE) {
  
  if (na.rm) {
    x <- stats::na.omit(x)
  }
  
  if (length(x) < 3) {
    stop("kurt() requires at least 3 observations.")
  }
  
  mu <- mean(x)
  xc <- x - mu
  
  s2 <- mean(xc^2)   # population variance
  
  if (s2 == 0) {
    stop("kurt() is undefined for zero-variance data.")
  }
  
  m4 <- mean(xc^4)   # fourth central moment
  
  k <- m4 / (s2^2)
  
  if (excess) {
    k <- k - 3
  }
  
  return(k)
}

#-------------------------------------------------------------------------------
# Compute octile skewness
#
# Robust skewness based on octiles (12.5%, 50%, 87.5%).
# Less sensitive to outliers than moment skewness.
# Useful for skewed soil and environmental data.
#-------------------------------------------------------------------------------
octile_skew <- function(x, na.rm = TRUE, type = 7) {
  
  if (na.rm) x <- stats::na.omit(x)
  
  n <- length(x)
  if (n < 3)
    stop("octile_skew() requires at least 3 observations.")
  
  q <- stats::quantile(x, probs = c(1/8, 1/2, 7/8),
                       names = FALSE, type = type)
  
  denom <- q[3] - q[1]
  
  if (denom == 0)
    return(0)  # no spread → symmetric by definition
  
  os <- ((q[3] - q[2]) - (q[2] - q[1])) / denom
  
  return(os)
}

#-------------------------------------------------------------------------------
# Internal: compute Tukey fences and outlier diagnostics
#-------------------------------------------------------------------------------
.tukey_stats <- function(x) {
  
  q1 <- stats::quantile(x, 0.25, names = FALSE)
  q2 <- stats::quantile(x, 0.50, names = FALSE)
  q3 <- stats::quantile(x, 0.75, names = FALSE)
  
  iqr <- q3 - q1
  
  inner_low  <- q1 - 1.5 * iqr
  inner_high <- q3 + 1.5 * iqr
  
  outer_low  <- q1 - 3 * iqr
  outer_high <- q3 + 3 * iqr
  
  list(
    q1 = q1, q2 = q2, q3 = q3,
    inner_low = inner_low,
    inner_high = inner_high,
    outer_low = outer_low,
    outer_high = outer_high,
    mild_outliers    = which(x < inner_low | x > inner_high),
    extreme_outliers = which(x < outer_low | x > outer_high)
  )
}

#-------------------------------------------------------------------------------
# Summary statistics for a numeric vector
#-------------------------------------------------------------------------------
summary_stats <- function(x, decimals = NULL, na.rm = TRUE, detailed_outliers = FALSE) {
  
  if (na.rm) x <- stats::na.omit(x)
  
  if (!is.numeric(x))
    stop("summary_stats() requires a numeric vector.")
  
  if (length(x) < 3)
    stop("summary_stats() requires at least 3 observations.")
  
  tk <- .tukey_stats(x)   # your internal Tukey stats function
  
  # ---------------------------
  # Core statistics
  # ---------------------------
  out <- data.frame(
    "Mean"            = mean(x),
    "Median"          = tk$q2,
    "Q1"              = tk$q1,
    "Q3"              = tk$q3,
    "Variance"        = stats::var(x),
    "SD"              = stats::sd(x),
    "Skewness"        = skew(x),
    "Octile Skewness" = octile_skew(x),
    "Kurtosis"        = kurt(x),
    "No. Outliers"    = length(union(tk$mild_outliers, tk$extreme_outliers)),
    check.names = FALSE
  )
  
  # ---------------------------
  # Optional detailed outliers
  # ---------------------------
  if (detailed_outliers) {
    out$N_Mild_Outliers    <- length(tk$mild_outliers)
    out$N_Extreme_Outliers <- length(tk$extreme_outliers)
    out$Inner_Lower_Fence  <- tk$inner_low
    out$Inner_Upper_Fence  <- tk$inner_high
    out$Outer_Lower_Fence  <- tk$outer_low
    out$Outer_Upper_Fence  <- tk$outer_high
  }
  
  # ---------------------------
  # Apply decimal places if requested
  # ---------------------------
  if (!is.null(decimals)) {
    out <- round(out, decimals)
  }
  
  return(out)
}

#-------------------------------------------------------------------------------
# Diagnostic plots:
#   Histogram + Tukey box overlay + Normal Q-Q plot
#-------------------------------------------------------------------------------
summary_plot <- function(x, varname = "x", na.rm = TRUE) {
  
  if (na.rm) x <- stats::na.omit(x)
  
  if (!is.numeric(x))
    stop("summary_plot requires a numeric vector.")
  
  tk <- .tukey_stats(x)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(1, 2))
  
  # ---------------------------
  # Histogram with Tukey box overlay
  # ---------------------------
  h <- hist(x, plot = FALSE)
  ymax <- max(h$counts)
  
  hist(x,
       col = "lightblue",
       border = "#000000",
       breaks = 30,
       main = "",
       xlab = varname,
       ylim = c(0, ymax * 1.25))
  
  box_ymin <- ymax * 1.10
  box_ymax <- ymax * 1.20
  box_mid  <- ymax * 1.15
  
  # box
  segments(tk$q1, box_ymin, tk$q3, box_ymin)
  segments(tk$q1, box_ymax, tk$q3, box_ymax)
  segments(tk$q1, box_ymin, tk$q1, box_ymax)
  segments(tk$q3, box_ymin, tk$q3, box_ymax)
  
  # whiskers (using inner fences)
  inliers <- x[x >= tk$inner_low & x <= tk$inner_high]
  lw <- min(inliers)
  uw <- max(inliers)
  
  segments(lw, box_mid, tk$q1, box_mid, lty = 2)
  segments(tk$q3, box_mid, uw, box_mid, lty = 2)
  
  # ---------------------------
  # Whisker end caps (adjacent values)
  # ---------------------------
  cap_half_height <- (box_ymax - box_ymin) * 0.15
  
  segments(lw, box_mid - cap_half_height, lw, box_mid + cap_half_height, lwd = 2.5)
  segments(uw, box_mid - cap_half_height, uw, box_mid + cap_half_height, lwd = 2.5)
  
  # median
  segments(tk$q2, box_ymin, tk$q2, box_ymax, lwd = 2)
  
  # ---------------------------
  # Mean (red thick)
  # ---------------------------
  mu <- mean(x)
  segments(mu, box_ymin, mu, box_ymax, lwd = 2, col = "red", lty = 2)
  
  
  # ---------------------------
  # Jittered raw data points
  # ---------------------------
  set.seed(1)  # reproducible jitter
  points(
    jitter(x, amount = diff(range(x)) * 0.01),
    runif(length(x), box_ymin, box_ymax),
    pch = 21,  # filled circle with border
    bg = adjustcolor("steelblue", alpha.f = 0.3),  # fill color with transparency
    col = "black",  # border color
    cex = .65
  )
  
  # inner fences (mild outliers)
  abline(v = tk$inner_low,  col = "orange", lty = 2)
  abline(v = tk$inner_high, col = "orange", lty = 2)
  
  # inner fences (extreme outliers)
  abline(v = tk$outer_low,  col = "red", lty = 2)
  abline(v = tk$outer_high, col = "red", lty = 2)
  
  # ---------------------------
  # Normal Q-Q plot
  # ---------------------------
  qq <- qqnorm(x, main = "", pch = 21, col = "black", bg = adjustcolor("steelblue", alpha.f = 0.3))
  qqline(x, col = "blue", lty = 2)
  
  # ---------------------------
  # Mild outliers (orange)
  # ---------------------------
  if (length(tk$mild_outliers) > 0) {
    points(qq$x[tk$mild_outliers],
           qq$y[tk$mild_outliers],
           pch = 21,
           bg = adjustcolor("orange", alpha.f = 0.8),
           col = "black")
  }
  
  # ---------------------------
  # Extreme outliers (red)
  # ---------------------------
  if (length(tk$extreme_outliers) > 0) {
    points(qq$x[tk$extreme_outliers],
           qq$y[tk$extreme_outliers],
           pch = 21,
           bg = adjustcolor("red", alpha.f = 0.8),
           col = "black")
  }
}

#-------------------------------------------------------------------------------
# Histogram with Tukey boxplot overlay and outlier fences
#
# Colors:
#   Orange → mild outlier limits (±1.5 IQR)
#   Red    → extreme outlier limits (±3 IQR)
#-------------------------------------------------------------------------------
hist_plot <- function(x, varname = "x", na.rm = TRUE) {
  
  # ---------------------------
  # Input checks
  # ---------------------------
  if (na.rm) x <- stats::na.omit(x)
  
  if (!is.numeric(x))
    stop("hist_plot() requires a numeric vector.")
  
  if (length(x) < 3)
    stop("hist_plot() requires at least 3 observations.")
  
  # ---------------------------
  # Tukey statistics
  # ---------------------------
  tk <- .tukey_stats(x)
  
  # ---------------------------
  # Histogram base
  # ---------------------------
  par(mfrow = c(1, 1))
  
  h <- hist(x, plot = FALSE)
  ymax <- max(h$counts)
  
  hist(x,
       col = "lightblue",
       border = "#000000",
       breaks = 30,
       main = "",
       xlab = varname,
       ylim = c(0, ymax * 1.25))
  
  # ---------------------------
  # Vertical span for boxplot inside histogram
  # ---------------------------
  box_ymin <- ymax * 1.10
  box_ymax <- ymax * 1.20
  box_mid  <- (box_ymin + box_ymax) / 2
  
  # ---------------------------
  # Draw IQR box
  # ---------------------------
  rect(tk$q1, box_ymin, tk$q3, box_ymax, border = "black", col = NA)
  
  # End caps at Q1 and Q3 (small horizontal ticks)
  segments(tk$q1, box_ymin, tk$q1, box_ymax, lwd = 2)
  segments(tk$q3, box_ymin, tk$q3, box_ymax, lwd = 2)
  
  # ---------------------------
  # Whiskers (non-outlier range)
  # ---------------------------
  inliers <- x[x >= tk$inner_low & x <= tk$inner_high]
  lw <- min(inliers)
  uw <- max(inliers)
  
  segments(lw, box_mid, tk$q1, box_mid, lwd = 1.5, lty = 2)
  segments(tk$q3, box_mid, uw, box_mid, lwd = 1.5, lty = 2)
  
  # ---------------------------
  # Whisker end caps (adjacent values)
  # ---------------------------
  cap_half_height <- (box_ymax - box_ymin) * 0.15
  
  segments(lw, box_mid - cap_half_height, lw, box_mid + cap_half_height, lwd = 2.5)
  segments(uw, box_mid - cap_half_height, uw, box_mid + cap_half_height, lwd = 2.5)
  
  # ---------------------------
  # Median (black thick)
  # ---------------------------
  segments(tk$q2, box_ymin, tk$q2, box_ymax, lwd = 2.5)
  
  # ---------------------------
  # Median (black thick)
  # ---------------------------
  segments(tk$q2, box_ymin, tk$q2, box_ymax, lwd = 2.5, col = "black")
  
  # ---------------------------
  # Mean (red thick)
  # ---------------------------
  mu <- mean(x)
  segments(mu, box_ymin, mu, box_ymax, lwd = 2, col = "red", lty = 2)
  
  
  # ---------------------------
  # Jittered raw data points
  # ---------------------------
  set.seed(1)  # reproducible jitter
  points(
    jitter(x, amount = diff(range(x)) * 0.01),
    runif(length(x), box_ymin, box_ymax),
    pch = 21,  # filled circle with border
    bg = adjustcolor("steelblue", alpha.f = 0.3),  # fill color with transparency
    col = "black",  # border color
    cex = 1
  )
  
  # ---------------------------
  # Mild outlier fences (orange)
  # ---------------------------
  abline(v = tk$inner_low,  col = "orange", lty = 2, lwd = 2)
  abline(v = tk$inner_high, col = "orange", lty = 2, lwd = 2)
  
  # ---------------------------
  # Extreme outlier fences (red)
  # ---------------------------
  abline(v = tk$outer_low,  col = "red", lty = 2, lwd = 2)
  abline(v = tk$outer_high, col = "red", lty = 2, lwd = 2)
}


################################################################################
histplot<-function(x,varname){

# plot a histogram with boxplot in x indicating
# any probable outliers by Tukey's criterion

x<-na.drop(x)
if(missing(varname)){varname<-"x"}

Q1<-quantile(x,prob=0.25)
Q3<-quantile(x,prob=0.75)
Q2<-quantile(x,prob=0.5)
hspread<-Q3-Q1
Fu<-Q3+3*hspread
Fl<-Q1-3*hspread


# ols,oll: values below and above outer fences
# posols,posoll: values below and above inner fences 
# (so ols and posols overlap, as do oll and posoll
#
ols<-which(x<Fl)
oll<-which(x>Fu)
posols<-which(x<(Q1-1.5*hspread))
if(length(posols)==0){
lw<-min(x)}else{
lw<-min(x[-posols])}
posoll<-which(x>(Q3+1.5*hspread))
if(length(posoll)==0){
uw<-max(x)}else{
uw<-max(x[-posoll])}

ol<-c(ols,oll) # combined outlier set

ymax<-max((hist(x,plot=F))$counts)
hist(x,main="",col="AliceBlue", xlab=varname,ylim=c(0,(ymax*1.25)))

boxmin<-ymax*1.1
boxmax<-ymax*1.2
boxmid<-ymax*1.15

lines(c(Q1,Q3),c(boxmin,boxmin))
lines(c(Q1,Q3),c(boxmax,boxmax))
lines(c(Q1,Q1),c(boxmin,boxmax))
lines(c(Q3,Q3),c(boxmin,boxmax))
lines(c(Q1,lw),c(boxmid,boxmid))
lines(c(Q3,uw),c(boxmid,boxmid))
lines(c(Q2,Q2),c(boxmin,boxmax),lwd=2)

lines(c(Fu,Fu),c(10,boxmid),lty=5,col="red")
lines(c(Fl,Fl),c(10,boxmid),lty=5,col="red")

}
###########################################################################

outliers<-function(x,trim){

# compute summary statistics of values in x

if(missing(trim)){trim<-F}

x<-na.drop(x)
Q1<-quantile(x,prob=0.25)
Q3<-quantile(x,prob=0.75)
Q2<-quantile(x,prob=0.5)
hspread<-Q3-Q1
Fu<-Q3+3*hspread
Fl<-Q1-3*hspread

# ols,oll: values below and above outer fences
# posols,posoll: values below and above inner fences 
# (so ols and posols overlap, as do oll and posoll
#

ols<-which(x<Fl)
oll<-which(x>Fu)
posols<-which(x<(Q1-1.5*hspread))
if(length(posols)==0){
lw<-min(x)}else{
lw<-min(x[-posols])}
posoll<-which(x>(Q3+1.5*hspread))
if(length(posoll)==0){
uw<-max(x)}else{
uw<-max(x[-posoll])}

ol<-c(ols,oll) # combined outlier set

nol<-length(ol)

print(paste(nol," value(s) lie outwith the outer fences (Tukey, 1977)"),quote=F)

if(nol!=0){
print(paste("Outlying values are:"),quote=F)
print(x[ol])
print(paste("Indices of outlying values are"),quote=F)

if(trim==F){
return(ol)}else{
print(ol)
return(x[-ol])
}
}
}


###############################################################
###############################################################

okgridvar<-function(space,modtype,c0,c1,a1){

x0<-c(0,0)
x<-rep(seq(-2.5,2.5,1),6)*space
y<-rep(seq(-2.5,2.5,1),rep(6,6))*space
xg<-cbind(x,y)

a<-matrix(1,37,37)
a[37,37]<-0
b<-rep(1,37)

for (i in 1:36){
lag<-ddis(i,xg)
b[i]<-sv(lag,modtype,c0,c1,a1)
for (j in i:36){
if(i==j){a[i,j]<-0.0
}else{
lag<-dis(i,j,xg)
a[i,j]<-sv(lag,modtype,c0,c1,a1)
a[j,i]<-a[i,j]}
}
}

ai<-solve(a)
lam<-ai%*%b
kv<-t(b)%*%lam
return(kv)
}

###############################################################
###############################################################
#
#  Functions below not yet documented
#
###############################################################
###############################################################


sv<-function(lag,modtype,c0,c1,a1){

hovera<-lag/a1

#modtype<-"Sph" or "Exp" 

if(modtype=="Exp"){semiv<-c0+c1*(1-exp(-hovera))
} else { 
if(lag>a1){semiv<-c0+c1}
if(lag<=a1){sph<-1.5*hovera-0.5*hovera^3
semiv<-c0+c1*sph}
}
return(semiv)
}

###############################################################
###############################################################

dis<-function(i,j,xg){

dx<-xg[i,1]-xg[j,1]
dy<-xg[i,2]-xg[j,2]
return(sqrt((dx*dx)+(dy*dy)))
}

ddis<-function(i,xg){

dx<-xg[i,1]
dy<-xg[i,2]
return(sqrt((dx*dx)+(dy*dy)))
}

ocskew<-function(x){

percs<-quantile(x,c(1/8,1/2,7/8))
o7<-percs[3]
o4<-percs[2]
o1<-percs[1]

os<-((o7-o4)-(o4-o1))/(o7-o1)
return(as.numeric(os))
}


#####################################################################
na.drop<-function(xin){
noNA<-as.numeric(length(which(is.na(xin)==T)))
if(noNA>0){
x<-as.numeric(na.omit(xin))
print(paste(noNA," missing value(s) removed"),quote=F)
}else{
x<-xin
}
return(x)
}

######################################################################


vanG<-function(ps1,id,xidep){
h<-xidep[,1]
thr<-ps1[id,1]
ths<-ps1[id,2]
alp<-ps1[id,3]
nscal<-ps1[id,4]
mscal<-ps1[id,5]
den<-(1+((abs(alp*h))^nscal))^mscal
return(thr+((ths-thr)/den))
}

vanGdraw<-function(psi,h){
thr<-psi[1]
ths<-psi[2]
alp<-psi[3]
nscal<-psi[4]
mscal<-psi[5]
den<-(1+((abs(alp*h))^nscal))^mscal
return(thr+((ths-thr)/den))
}

campB<-function(ps1,id,xidep){
h<-xidep[,1]
ths<-ps1[id,1]
alp<-ps1[id,2]
nscal<-ps1[id,3]
return(ths*((alp*h)^(-1/nscal)))
}


campBdraw<-function(psi,h){

ths<-psi[1]
alp<-psi[2]
nscal<-psi[3]

return(ths*((alp*h)^(-1/nscal)))
}



plotgroups<-function(x,y,icl,data.df,xlab,ylab,xylog){
# function to plot variables x and y in dataframe data.df
# groups with different variables (in icl) are plotted with contrasting
# symbols.  xlab and ylab are optional labels.  Set xylog to x or
# y or xy for log scales, n for neither

options(scipen=999)
if(missing(xlab)){xlab<-"x"}
if(missing(ylab)){ylab<-"y"}
if(missing(xylog)){xylog<-"n"}

cnames<-colnames(data.df)
indl<-which(cnames==icl)
indx<-which(cnames==x)
indy<-which(cnames==y)

levels(data.df[,indx])<-seq(1,nlevels(data.df[,indx]))


if(xylog=="n"){
plot(data.df[,indx],data.df[,indy],xlab=xlab,ylab=ylab,pch=16,col=data.df[,indl])
}else{
plot(data.df[,indx],data.df[,indy],log=xylog,xlab=xlab,ylab=ylab,pch=16,col=data.df[,indl])
}

}


mean_censor<-function(y,cen,log.t){
if(missing(log.t)){log.t=F}

mean.guess<-mean(y)
sd.guess<-sd(y)

ncen<-length(which(is.na(censor(y,cen))))
oo<-optim(c(mean.guess,sd.guess),nllcen,cen=cen,y=y)
oo2<-optim(oo$par,nllcen,cen=cen,y=y,hessian=T)

se<-sqrt((solve(oo2$hessian))[1,1])

mean.est<-oo2$par[1]
sd.est<-oo2$par[2]

# approximate df by n-1

nobs<-length(y)
tval<-qt(0.975,(nobs-1))
clu<-oo2$par[1]+tval*se
cll<-oo2$par[1]-tval*se

mean.naive<-mean(censor(y,cen),na.rm=T)

op1<-matrix(c((-1*oo2$value),ncen,mean.est,sd.est,cll,clu,mean.naive),1,7)
colnames(op1)<-c("log_likelihood","Number_censored","Mean","SD","Lower_95%","Upper_95%","Mean_>_DL_only")



if(log.t==T){

median.unbiased<-exp(mean.est)
cll.bt<-exp(cll)
clu.bt<-exp(clu)


op2<-matrix(c(median.unbiased,cll.bt,clu.bt),1,3)
colnames(op2)<-c("Back_median","Back_Lower_95%","Back_Upper_95%")
op1<-cbind(op1,op2)
}

return(op1)
}

#
censor<-function(x,cen){
censored<-which(x<cen)
x[censored]<-NA
return(x)
}

nllcen<-function(theta,cen,y){
#
#
#
mu<-theta[1]
sig<-theta[2]

y<-censor(y,cen)
n_c<-length(which(is.na(y)))
Phi_c<-log(pnorm(cen,mu,sig))
z<-y[-which(is.na(y))]

dens<-dnorm(z,mu,sig,log=T)

return(-n_c*Phi_c-sum(dens))
}



