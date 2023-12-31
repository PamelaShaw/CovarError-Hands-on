library(MASS)
library(boot)
library(mvtnorm)

library(devtools)
#> Loading required package: usethis
#install_github("lboe23/sandwich2stage", subdir="pkg")
library(sandwich2stage)

### Set up some simulated data and fix seed for reproducibility
	set.seed(353)
	n <- 1000  ### Cohort sample size
	n1 <- 100  ### Calibration subset sample size

### anti-logit utility function
	expit=function(x){exp(x)/(1+exp(x))}
	
	xzcorr <- 0.25
	XZ <- rmvnorm(n, mean=c(0,0), sigma=array(c(1,xzcorr,xzcorr,1),dim=c(2,2)))
	x <- XZ[,1]
	z <- XZ[,2]

	y <- rbinom(n,1,expit(-1+log(1.5)*x+log(1.75)*z))	  ###logistic regression model
	mean(y)

	reliability <- 0.5	### Attenuation coef
	sigma_u_sq <- 1/reliability - 1
	xstar1 <- x+rnorm(n, sd=sigma_u_sq^0.5)
	xstar2 <- x+rnorm(n, sd=sigma_u_sq^0.5)

	mydat <- data.frame(y,z,xstar1,xstar2)

	###The gold standard analysis and the naive analysis

	#naive model
	##Suppose have xstar1,xstar2 on everyone, where xstar=X+u
	mydat$xstarbar<-(mydat$xstar1+mydat$xstar2)/2
	naive.model=glm(y~xstarbar,family="binomial",data=mydat)

	summary(naive.model)

	#### True Model
	true.model=glm(y~x,family="binomial",data=mydat)
	summary(true.model)
	
###Regression calibration model

  varU<-var(mydat$xstar1-mydat$xstar2)/2
  meanXhat<-mean(mydat$xstarbar)
  lambda<-(var(mydat$xstarbar) -varU/2)/ var(mydat$xstarbar)
  mydat$xhat<- lambda*mydat$xstarbar + (1-lambda)*meanXhat

#perform regression calibration
final.model=glm(y~xhat,family="binomial",data=mydat)
summary(final.model)$coef

###Need to fix Standard errors for final model

bootstrap.function<-function(dat,inds){
	mydat.boot<-dat[inds,]
	varU<-var(mydat.boot$xstar1-mydat.boot$xstar2)/2
	meanXhat<-mean(mydat.boot$xstarbar)
	lambda<-(var(mydat.boot$xstarbar) -varU/2)/ var(mydat.boot$xstarbar)
	mydat.boot$xhat<- lambda*mydat.boot$xstarbar + (1-lambda)*meanXhat
	
	final.model=glm(y~xhat,family="binomial",data=mydat.boot)
	return(final.model$coef)
}
my.boot<-boot(mydat,bootstrap.function,R=1000)
my.boot

	
	##### Alternative approach that makes it easy to incorporate covariates, easier but perhaps less efficient.
	##### Works if you have repeat meausres or just a classical measurement error in the subset
	###Regression calibration model
	#perform regression calibration
	rcfit<-lm(xstar2~xstar1+z,data=mydat)
	mydat$xhat<-predict(rcfit,newdata=mydat)
	final.model=glm(y~xhat+z,family="binomial",data=mydat)

	summary(final.model)
	
###Standard errors

bootstrap.function<-function(dat,inds){
	mydat.boot<-dat[inds,]
	rcfit.boot<-lm(xstar2~xstar1+z,data=mydat.boot)
	mydat.boot$xhat<-predict(rcfit.boot,newdata=mydat.boot)
	
	final.model=glm(y~xhat + z,family="binomial",data=mydat.boot)
	return(final.model$coef)
}
my.boot<-boot(mydat,bootstrap.function,R=1000)
bsSD<- apply(my.boot$t,2,sd)
bsSD

t.stat<-coef(final.model)/bsSD
t.stat

### Can calculate sandwich SE for comparison

### Need to use survey package for sandwich2stage

### Declare a simple random sample design for svyglm()
names(mydata)
datadesign <- svydesign(id=~1, data=mydat)

## calibration model is the stage 1 model
stage1.model<-svyglm(xstar2~xstar1+z,design=datadesign,family=gaussian(),data=mydat)
datadesign <- update(datadesign,xhat =predict(stage1.model,newdata=datadesign$variables),ID=1:nrow(mydat))

### outcome model is the stage 2 model
stage2.model<-  svyglm(y ~ xhat+z,design=datadesign,family=binomial())
sandwich.object<-sandwich2stage(stage1.model,stage2.model,xstar="xstar1",xhat="xhat",Stage1ID="ID",Stage2ID="ID")
swvar<-vcov(sandwich.object) 
### swvar is 6 by 6 var-covar matrix (3 calibration params and 3 outcome model params )
swvar 
### These are the SE for the three parameters in the stage2.model

sqrt(diag(swvar)[4:6])




