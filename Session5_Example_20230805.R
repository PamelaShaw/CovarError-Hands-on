library(openxlsx)
library(survival)
library(simexaft)  ### simex for the accelerated failure time model, PH survival
library(boot)

#> Loading required package: use this to install 
#install_github("lboe23/sandwich2stage", subdir="pkg")
library(sandwich2stage)


### Inspired by WHI dietary data 
data<-read.xlsx("~/Library/CloudStorage/OneDrive-KaiserPermanente/Meetings/JSM 2023/simulated_data.xlsx")
names(data)

### true log energy intake = logEn
### selfreported log energy intake = logEn_sr (systematic measurement error)
### biomarker log energy intake = logEn_bm  (classical measurement error)
### repeat measures of these variables have a 2 after them
### Physical activity = step_1000 (units 1000 steps)
### Body mass index = BMI
### age at baseline = age
### race = race


### Note sandwich2stage and simexFit did not seem to like factor variables, so create the binary dummary variables for each race
data$black<-ifelse(data$race=="black",1,0)
data$other<-ifelse(data$race=="other",1,0)
data$white<-ifelse(data$race=="white",1,0)

### turn race into a factor and make the biggest racial group the reference category
data$myrace<-relevel(factor(data$race),ref="white")
true<-coxph(Surv(time,delta)~logEn+age+BMI+step_1000+myrace, data=data)
summary(true)

### Fit ignoring meas error
naive<-coxph(Surv(time,delta)~logEn_sr+age+BMI+step_1000+myrace, data=data)
summary(naive)

### Fit calibration
calfit<-lm(logEn_bm~logEn_sr + age+BMI+step_1000+ myrace,data=data)
data$logEn.hat<- predict(calfit,newdata=data)
summary(calfit)

### Perform regression calibration
rc<-coxph(Surv(time,delta)~logEn.hat+age+BMI+step_1000+myrace, data=data)
summary(rc)

### Get bootstrapped SE and confidence intervals
set.seed(806)
bootstrap.function<-function(dat,inds){
    mydat.boot<-dat[inds,]
    rc.model=lm(logEn_bm~ logEn_sr+ age + BMI+step_1000+myrace,data=mydat.boot)
    xhat=predict(rc.model,newdata=mydat.boot)
    final.model=coxph(Surv(time,delta)~logEn.hat+age+BMI+step_1000+myrace, data=mydat.boot)
    return(final.model$coef)
}
data$subset<-ifelse(data$subset==1,1,0)  ### Strata variable needed for bootsrap
my.boot<-boot(data,bootstrap.function, strata=data$subset,R=250) ### R low to speed up for class, I would choose R=500 or 1000
my.boot
boot.ci(my.boot,type=c("norm","perc"))

### can get sandwich SE
### Declare a simple random sample design for svyglm()
names(data)
datadesign <- svydesign(id=~1, data=data)

## calibration model is the stage 1 model
stage1.model<-svyglm(logEn_bm~logEn_sr+age+BMI+step_1000+black+other, design=datadesign, family=gaussian(),data=data)
datadesign <- update(datadesign, xhat =predict(stage1.model,newdata=datadesign$variables),ID=1:nrow(data))

### outcome model is the stage 2 model
stage2.model<-  svycoxph(Surv(time,delta)~ xhat+age+BMI+step_1000+black+other, design=datadesign)

sandwich.object<-sandwich2stage(stage1.model,stage2.model,xstar="logEn_sr",xhat="xhat",Stage1ID="ID",Stage2ID="ID")
swvar<-vcov(sandwich.object) 
### swvar is 13 by 13 var-covar matrix (3 calibration params and 3 outcome model params )
swvar 
### These are the SE for the 6 parameters in the stage2.model
sqrt(diag(swvar)[8:13])



##############################Explore SIMEX#############################
#### Simex assumes classical measruement error and a known error function
#### have the systematic error variable on everyone. 
#### Have the classical measurment error variable on subset, but very few events on that subset

###fit a AFT model and use simex for accelerated failure time models-- simexaft package
set.seed(120) 
### Dont have repeats of the self reported diet, so could guess the measurement error variance is = to total variance in exposure
ind <- c("logEn_sr") 
err.mat <- var(data$logEn_sr)


### saving some typing by creating formula variables
true.formula<- "Surv(time, delta) ~ logEn_sr + age + BMI + step_1000+black+other"
naive.formula<- "Surv(time, delta) ~ logEn_sr + age + BMI + step_1000+black+other"

## Since data are comptuer generated we can consider the true exposure, for a benchmark
 trueAFT<-survreg(formula = formula(true.formula), data = data, dist = "weibull")

#### Naive model using errorprone self-reported energy with no correction for error
naiveAFT<-survreg(formula = formula(naive.formula), data = data, dist = "weibull")

### Look at SIMEX for weibull model, using quadratic exptrapolation
simexFit <- simexaft(formula = formula(naive.formula), data = data, SIMEXvariable = ind, repeated = FALSE, err.mat = err.mat, B = 50, lambda = seq(0, 2, 0.1),extrapolation = "quadratic", dist = "weibull") 
summary(simexFit)
plotsimexaft(simexFit,"logEn_sr",ylimit=c(-3,1),extrapolation="quadratic")

### Try linear extrapolation
simexFitLin <- simexaft(formula = formula(naiveAFT), data = data, SIMEXvariable = ind, repeated = FALSE, err.mat = err.mat, B = 50, lambda = seq(0, 2, 0.1),extrapolation = "linear", dist = "weibull") 
summary(simexFitLin)
plotsimexaft(simexFitLin,"logEn_sr",ylimit=c(-3,1),extrapolation="linear")

### Hazard ratio (HR) estimated by true Cox model
trueCox<-coxph(as.formula(true.formula), data=data)
exp(trueCox$coef)

###  HR estimated by true AFT model, similar to Cox model
exp(-trueAFT$coef/trueAFT$scale)

### HR estimated by simex AFT model: Does terribly because systematic and not classical error in the self-reported energy variable
### Some inflation here, due to the error being more complicated than the simple classical error SIMEX assumes. Also had to guess at error variance.
exp(-simexFit$coef/simexFit$scale)

##########
#### Look at sensitiivty to the assumed error variance.
##########
err.mat2<-.5*err.mat
simexFit2 <- simexaft(formula = formula(naive.formula), data = data, SIMEXvariable = ind, repeated = FALSE, err.mat = err.mat2, B = 50, lambda = seq(0, 2, 0.1),extrapolation = "quadratic", dist = "weibull") 
exp(-simexFit2$coef/simexFit2$scale)

##### Could look at actual error variance, but simex error model still wrong.
err.mat3<- var(data$logEn - data$logEn_sr)
simexFit3<- simexaft(formula = formula(naive.formula), data = data, SIMEXvariable = ind, repeated = FALSE, err.mat = err.mat3, B = 50, lambda = seq(0, 2, 0.1),extrapolation = "quadratic", dist = "weibull") 
exp(-simexFit3$coef/simexFit3$scale)



