source("dplfit.R")

library(ggplot2)

#
# Helpful functions
#

# make_discrete_cdf - create data frame with points for the top and bottom of
# discontinuities in a discrete cumulative distribution function
make_discrete_cdf = function(xs, ys) {
  return(data.frame(
    x=rep(xs,each=2)[1:(2*length(xs)-1)],
    y=c(1,rep(ys[1:(length(ys)-1)],each=2)))) }

# make_ecdf_df - create a dataframe in the same form as make_discrete_cdf for
# the cdf of an empirical sample
make_ecdf_df = function(vs) {
  ecd = ecdf(vs)
  xs = environment(ecd)$x
  ys = environment(ecd)$y
  return(make_discrete_cdf(xs, 1-ys)) }


# Tests
within_se = function (val, sedf) {
  return(with(sedf, val > alpha - 1.96*se & val < alpha + 1.96*se)) }

n_replicates = 1000
estimate_pareto = table(replicate(n_replicates,
  within_se(2,paretofit_se(rpareto(1000,2,1),1))))
pareto_empirical_ci = estimate_pareto[2]/n_replicates
pareto_empirical_ci
estimate_pareto = table(replicate(n_replicates,
  within_se(2,dplfit_bin_se(rpareto(1000,2,1),2,1))))
pareto_empirical_ci = estimate_pareto[2]/n_replicates
pareto_empirical_ci
estimate_pareto = table(replicate(n_replicates,
  within_se(2,dplfit_bin_se(rpareto(1000,2,1),2,1))))
pareto_empirical_ci = estimate_pareto[2]/n_replicates
pareto_empirical_ci

#
# Figure 1 - Tail distributions of powere laws
#
xs = 20^(seq(0,100)/100)
zeta_xs = 1:20
cdfs = rbind(
  data.frame(src="Pareto (continuous)", x=xs, y=ppareto(xs,2)),
  cbind(src="Zeta (discrete)", make_discrete_cdf(zeta_xs,pzeta(zeta_xs,2))),
  cbind(src="self-similar discrete", make_discrete_cdf(xs,pdpl(xs,2))))

ggplot(cdfs,aes(x=x,y=1-y)) + geom_line() + 
  coord_cartesian(expand=FALSE) + scale_y_log10() + scale_x_log10() +
  labs(y="log prob. observ. exceeds x", x="log x") + 
  facet_wrap(~src,ncol=1,scales="fixed")


#
# Figure 2a - Pareto and discrete power law fits.
#
xm = 1
lam = sqrt(2)
alpha = 2

set.seed(2)

ss <- rdpl(10000,alpha,lam,xm)
fits = function (xms) { 
  fits = data.frame()
  for(i in 1:length(xms)) {
    fits = rbind(fits, paretofit_se(ss,xms[i])) }
  return(cbind(xm=xms, fits)) }
ssdfits = function (xms) { 
  fits = data.frame()
  for(i in 1:length(xms)) {
    fits = rbind(fits, dplfit_se_(ss,lam,xms[i])) }
  return(cbind(xm=xms, fits)) }

xms <- xm * lam^seq(0,7)
indata = fits(xms)

nsegs = 40 # tunable graph-drawing precision
outdata = expand.grid(
  xm = xms[2:length(xms)],
  ind = seq(0,nsegs-1))
outdata$odxm = with(outdata, xm * lam^(-ind/nsegs))
outdata = cbind(outdata,fits(outdata$odxm)[,c("alpha","se")])

fitdata = ssdfits(xms)

xmlambdalogscale = 
  scale_x_log10(breaks=xms,
    labels=c(expression(x[m]), expression(x[m]* lambda),
      expression(x[m]* lambda^2), expression(x[m]* lambda^3),
      expression(x[m]* lambda^4), expression(x[m]* lambda^5),
      expression(x[m]* lambda^6), expression(x[m]* lambda^7)))

ggplot(indata,aes(x=xm,y=alpha,ymin=alpha-1.96*se,ymax=alpha+1.96*se)) + 
  labs(y=expression(paste("Estimated ", hat(alpha))),
    x=expression(paste("Assumed ",x[m]))) +
  # True alpha:
  geom_hline(yintercept=2) +
  geom_point() +
  geom_errorbar(width=0.02) + 
  geom_point(data=fitdata,shape=5) + 
  geom_errorbar(data=fitdata, width=0.02) + 
  geom_line(data=outdata,aes(x=odxm,y=alpha,group=xm)) +
  geom_line(data=outdata,aes(x=odxm,y=(alpha+1.96*se),group=xm),linetype="11")+
  geom_line(data=outdata,aes(x=odxm,y=(alpha-1.96*se),group=xm),linetype="11")+
  geom_ribbon(data=outdata,
    aes(x=odxm,ymin=(alpha-1.96*se),ymax=(alpha+1.96*se),group=xm),
    fill="#000000", alpha=0.1,) +
  scale_y_continuous(breaks=c(0,1,2,3,4),
    lim=c(0,with(outdata,max(alpha+1.96*se)))) +
  xmlambdalogscale

#
# Figure 2b - Imperfect branching sample empirical cdf.
#

# Generate branching tree with OCaml code
system("ocamlbuild -pkgs gsl data.native")
ss <- as.numeric(system(intern=TRUE,
  "./data.native -t randombeta -a 2. -i 90.509668 -e 2.5 -x 1."))

# Randomly sample with replacement from the synthetic tree
ss = sample(ss, 10000) 

# Plot as in Fig 1.
cdfs = cbind(src="imperfect branching", make_ecdf_df(ss))
ggplot(cdfs,aes(x=x,y=y)) + geom_line() + 
  coord_cartesian(expand=FALSE) + scale_y_log10() + scale_x_log10() +
  labs(y=expression("log pr." >= x), x="log x") + 
  facet_wrap(~src,ncol=1,scales="fixed")

#
# Figure 2b - Pareto and discrete power law fits.
#
xm = 1
lam = sqrt(2)
alpha = 2

set.seed(2)

ss <- as.numeric(system(intern=TRUE,
  "./data.native -t randombeta -a 2. -i 100 -e 2.5 -x 1."))
ss = sample(ss, 10000)
fits = function (xms) { 
  fits = data.frame()
  for(i in 1:length(xms)) {
    fits = rbind(fits, paretofit_se(ss,xms[i])) }
  return(cbind(xm=xms, fits)) }
ssdfits = function (xms) {
  fits = data.frame()
  for(i in 1:length(xms)) {
    fits = rbind(fits, dplfit_bin_se(ss,lam,xms[i])) }
  return(cbind(xm=xms, fits)) }

xms <- xm * lam^seq(0,7)
indata = fits(xms)

nsegs = 40 # tunable graph-drawing precision
outdata = expand.grid(
  xm = xms[length(xms)],
  ind = seq(0,7*nsegs))
outdata$odxm = with(outdata, xm * lam^(-ind/nsegs))
outdata = cbind(outdata,fits(outdata$odxm)[,c("alpha","se")])

fitdata = ssdfits(xms)

ggplot(indata,aes(x=xm,y=alpha,ymin=alpha-1.96*se,ymax=alpha+1.96*se)) + 
  labs(y=expression(paste("Estimated ", hat(alpha))),
    x=expression(paste("Assumed ",x[m]))) +
  geom_hline(yintercept=2) +
  geom_point(data=fitdata,shape=5) + 
  geom_errorbar(data=fitdata, width=0.04) + 
  geom_line(data=outdata,aes(x=odxm,y=alpha,group=xm)) +
  geom_line(data=outdata,aes(x=odxm,y=(alpha+1.96*se),group=xm),linetype="11")+
  geom_line(data=outdata,aes(x=odxm,y=(alpha-1.96*se),group=xm),linetype="11")+
  geom_ribbon(data=outdata,
    aes(x=odxm,ymin=(alpha-1.96*se),ymax=(alpha+1.96*se),group=xm),
    fill="#000000", alpha=0.1,) +
  scale_y_continuous(breaks=c(1.5,2.0,2.5),
    lim=with(outdata,c(min(alpha-1.96*se),max(alpha+1.96*se)))) +
  xmlambdalogscale

#
# Figure 3. Earthquake and Bronchia
#
# 3a. Earthquake figure:
#
mags <- read.delim("data/earthcat/mags.all",header=FALSE)$V1
mags = mags[mags > 0]
mags = 10^mags
if(length(mags) < 1) {
  warning("Earthquake mags.all missing. Please run earthcat/make_earthcat.sh") }

cdf_lm_fit_res = function(xs,xm) {
  ecd = ecdf(xs[xs >= xm & !is.na(xs)])
  xs = environment(ecd)$x
  ys = environment(ecd)$y
  res = summary(lm(log(1-ys[1:(length(ys)-1)]) 
    ~ log(xs[1:(length(xs)-1)])))
  return(res) }

paretoslope = paretofit(mags, 10^2)
ssdslope = dplfit_(mags, 10^0.1, 10^2)
magscdf = cdf_lm_fit_res(mags, 10^2)

xmoffset = log(1-length(mags[mags < 10^2 & !is.na(mags)])/length(mags))
basemult = log(2,10)/log(2)

fits = rbind(
  data.frame(est = "Continuous MLE",
    intercept=paretoslope*log(10^2,10)+xmoffset*basemult,
    slope=-paretoslope),
  data.frame(est = "Discrete MLE",
    intercept=ssdslope*log(10^2,10)+xmoffset*basemult,
    slope=-ssdslope),
  data.frame(est = "Regression to emprical distribution",
    intercept=(magscdf$coefficients[1]+log(10^0.1) + xmoffset)*basemult,
    slope=magscdf$coefficients[2]))

eqw = 70
ggplot(make_ecdf_df(mags),aes(x=log(x,10),y=log(y,10))) + 
  coord_fixed(ratio=1,xlim=c(1.2,7.2)) +
  geom_abline(data=fits,
  aes(intercept=intercept,slope=slope,color=est)) +
  geom_line() + 
  scale_x_continuous(breaks=c(2,3,4,5,6,7,8),
    labels=c("2.0","3.0","4.0","5.0","6.0","7.0","8.0")) + 
  labs(y="log Prob(obs. > m)", x="Richter Magnitude, m")

#
# Figure 3b - Bronchial data
#

# Diameters < 2.0 mm are a 10% sample.
yehD6 <- read.delim("data/yeh1976tracheobronchial.D6.tsv", header=FALSE)
names(yehD6) <- c("diam")
yehD6$record = seq(1,1097)
yehD6 = yehD6[with(yehD6, c(yehD6[diam >= 20,"record"],
  rep(yehD6[diam < 20,"record"], 10))),]
yehD6 = yehD6[!is.na(yehD6$diam),]

yehlam = 2
yehxmin = 16

paretoslope = paretofit(yehD6$diam, yehxmin)
ssdslope = dplfit_bin(yehD6$diam, yehlam, yehxmin)
ssddata = log_bin(yehD6$diam, yehlam, yehxmin)
yehcdf = cdf_lm_fit_res(yehD6$diam, yehxmin)

xmoffset = log(1-length(yehD6$diam[yehD6$diam < yehxmin & !is.na(yehD6$diam)])/length(yehD6$diam))
basemult = log(2,10)/log(2)

fits = rbind(
  data.frame(est = "Continuous MLE",
    intercept=paretoslope*log(yehxmin,10),
    slope=-paretoslope),
  data.frame(est = "Discrete MLE",
    intercept=ssdslope*log(yehxmin,10),
    slope=-ssdslope),
  data.frame(est = "Regression to emprical distribution",
    intercept=(yehcdf$coefficients[1])*basemult,
    slope=yehcdf$coefficients[2]))

ggplot(rbind(data.frame(x=0,y=1),make_ecdf_df(yehD6$diam[yehD6$diam>=yehxmin])),
    aes(x=log(x,10),y=log(y,10))) + 
  coord_fixed(ratio=1) + 
  scale_x_continuous(breaks=log(c(20,50,100),10),labels=c(2,5,10))  +
  #geom_line(data=make_ssd_cdf(ssdslope,yehlam,yehxmin,yehnorders),alpha=0.5) +
  geom_abline(data=fits,
    aes(intercept=intercept,slope=slope,color=est)) +
  geom_line() + 
  geom_line(data=make_ecdf_df(ssddata),linetype="11") +
  labs(y="log Prob(obs. > m)", x="Diameter (mm)")
