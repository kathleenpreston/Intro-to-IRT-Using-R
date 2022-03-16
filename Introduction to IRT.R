#install.packages("psych")
#install.packages("GPArotation")
#install.packages("sirt")
#install.packages("ltm")
#install.packages("irtoys")
#install.packages("mirt")
#install.packages("latticeExtra")

  
library(psych)
library(GPArotation)
library(sirt)
library(msm)
library(ltm)
library(irtoys)
library(mirt)
library(latticeExtra)

# Read in Raw Data #
data<- read.delim("C:/Users/kpreston/Dropbox (Personal)/Research Projects/Talks/WPA Stats Talk/2018/Self-esteem Dichotomous.dat", header=FALSE)
data[data==9] <-NA #Replace missing code 9 with NA
nitems<-ncol(data)

#recode data
keys <- c(1,-1,1,1,-1,-1,1,-1,-1,1)
data <-reverse.code(keys,data,mini=rep(0,nitems),maxi=rep(1,nitems))

# Check response frequencies
response.frequencies(data)

(rfac <- tetrachoric(data,na.rm=TRUE))

## Assessing Unidimensionality via ECV, Reise (2012)
rfac <- tetrachoric(data,na.rm=TRUE)$rho
(g3 <- schmid(rfac, nfactors = 3, fm = "minres",digits=3,rotate="oblimin"))
(sumout <- colSums(g3$sl[,1:(ncol(g3$sl)-3)])^2)
(ECV2 <- sumout[1] / sum(sumout))

# Perform the Rasch model with ltm(), and save the results in modRasch
modRasch <- rasch(data, IRT.param = TRUE,constraint=cbind(ncol(data)+1,1))

round(coef(modRasch),2)  # Obtain difficulty and discrimination parameter estimates

plot(modRasch, type = "ICC", lwd = 2,legend=TRUE) #Item Characteristic Curves
plot(modRasch, type = "ICC", lwd = 2,legend=TRUE,item=3) #Only Item 3 ICC
plot(trf(est(data, model = "1PL",rasch="T", engine = "ltm"))) #Test Response Curve
plot(modRasch, type = "IIC", lwd = 2,legend = TRUE) #Item Information Curve
plot(modRasch, type = "IIC", lwd = 2,legend = TRUE,item=0) #Test Information Curve

summary(modRasch)

thetaRasch <-ability(data,est(data,model="1PL",rasch="T",engine="ltm"),method="BME")
round(head(thetaRasch),2)

item.fit(modRasch)


# Perform the 1PL analysis with ltm(), and save the results in "mod1PL"
mod1PL <- rasch(data, IRT.param = TRUE)
theta1PL <-ability(data,est(data,model="1PL",engine="ltm"),method="BME")
round(head(theta1PL),2)
round(coef(mod1PL),2)  # Obtain difficulty and discrimination parameter estimates
summary(mod1PL)
plot(mod1PL, type = "ICC", lwd = 2,legend=TRUE) #Item Characteristic Curves
plot(trf(est(data, model = "1PL", engine = "ltm"))) #Test Response Curve
plot(mod1PL, type = "IIC", lwd = 2,legend = TRUE) #Item Information Curve
plot(mod1PL, type = "IIC", lwd = 2,legend = TRUE,item=0) #Test Information Curve

item.fit(mod1PL)
# Compare Model fit
anova.rasch(modRasch,mod1PL)


# Perform the 2PL analysis with ltm(), and save the results in "mod2PL"
mod2PL <- ltm(data ~ z1, IRT.param = TRUE)
round(coef(mod2PL),2)  # Obtain difficulty and discrimination parameter estimates
theta2PL <-ability(data,est(data,model="2PL",engine="ltm"),method="BME")
round(head(theta2PL),2)
summary(mod2PL)
plot(mod2PL, type = "ICC", lwd = 2,legend=TRUE) #Item Characteristic Curves
plot(trf(est(data, model = "2PL", engine = "ltm"))) #Test Response Curve
plot(mod2PL, type = "IIC", lwd = 2,legend = TRUE) #Item Information Curve
plot(mod2PL, type = "IIC", lwd = 2,legend = TRUE,item=0) #Test Information Curve

item.fit(mod2PL)
# Compare Model fit
anova.rasch(mod1PL,mod2PL)



# Perform the 3PL analysis with ltm(), and save the results in "mod3PL"
mod3PL <- tpm(data, type="latent.trait",IRT.param = TRUE,control=c(optimizer="nlminb"))
round(coef(mod3PL),2)  # Obtain difficulty and discrimination parameter estimates
plot(mod3PL, type = "ICC", lwd = 2,legend=TRUE) #Item Characteristic Curves
plot(mod3PL, type = "IIC", lwd = 2,legend = TRUE) #Item Information Curve
plot(mod3PL, type = "IIC", lwd = 2,legend = TRUE,item=0) #Test Information Curve


## Confirming unidimensionality
unidimTest(mod2PL)

## Assessing Local Dependence using Yen's Q3 statistic, (Yen, 1984)
mod.q3 <- Q3( dat = data, theta = theta2PL[,1] , b = coef(mod2PL)[,1])
mod.q3$q3.matrix[upper.tri(mod.q3$q3.matrix,diag=TRUE)] <-0
round(mod.q3$q3.matrix,2)

#Criteria relative to average observed residual correlation from Christensen, Maransky & Horton (2017)
mod.q3$Q3.stat
mod.q3$Q3.stat[9]-abs(mod.q3$Q3.stat[1])


## Computing theta
theta <- seq(-4,4,.001)
apar <- coef(mod2PL)[,2]
bpar <- coef(mod2PL)[,1]

P <- matrix(0,nitems,length(theta))
L <- matrix(0,nitems,length(theta))

dis <- dnorm(theta)
dis <- dis/sum(dis)

person1 <- data[1,]

for(m in 1:nitems){
  P[m,] <- 1/(1+exp(-(apar[m]*(theta-bpar[m]))))
}

for(k in 1:length(person1)){
  if(person1[k] == 1){ L[k,] <- P[k,] }
  else if(person1[k] == 0) {L[k,] <- 1-P[k,] }
}
likelihood <- dis
for(j in 1:length(person1)){
  likelihood <- likelihood*L[j,]
}

plot(theta,likelihood,type="l")
( bme <- theta[which.max(likelihood)])
round(theta2PL[1,],2)

sumScores <- rowSums(data)
cor(theta2PL[,1],sumScores,use="pairwise.complete.obs")

