setwd('/Users/Erin/Documents/Classes/826/hw2/')
#####################################################
# From PROGRAM 12.1 in H&R
# Descriptive statistics from NHEFS data (Table 12.1)
#####################################################
 
# reading external data file
#install.packages("sas7bdat")
library(sas7bdat)
nhefs <- read.sas7bdat("data/nhefs_book.sas7bdat")
nhefs <- as.data.frame(nhefs)
head(nhefs)

# Preliminary testing indicated that sex, age, and race are good
# candidates for 
nhefs.use <- with(nhefs,data.frame(seqn,qsmk,wt82_71,sex,age,race,income,
                                   marital,school,smokeyrs))
nhefs.measures <- nhefs.use[,2:10]
nhefs.cor <- cor(nhefs.measures,use="pairwise.complete.obs",method="spearman")
cor.test(nhefs.cor)
with(nhefs.final,summary(qsmk))


# (a) Use a logistic regression to estimate a propensity score
#     for quitting smoking. Examine if there are some overlap 
#     between the propensity scores of the two groups (smoking 
#     and non-smoking) based on the model. Include the summary 
#     for your logistic regression and the overlap plot for your 
#     logistic regression.


# Trying out a few
glm1 <- glm(qsmk~
              age+factor(sex)+factor(race)+factor(income)+factor(marital)+school,
            data=nhefs.use,family="binomial")
drop1(glm1,test="LRT")
glm2 <- glm(qsmk~
              age+factor(sex)+factor(race)+factor(income)+school,
            data=nhefs.use,family="binomial")
drop1(glm2,test="LRT")
glm3 <- glm(qsmk~
              age+factor(sex)+factor(race)+factor(income),
            data=nhefs.use,family="binomial")
drop1(glm3,test="LRT")
glm4 <- glm(qsmk~
              age+factor(sex)+factor(race),
            data=nhefs.use,family="binomial")
drop1(glm4,test="LRT")

nhefs.final <- with(nhefs.use,
                    na.omit(data.frame(seqn,qsmk,age,sex,race,wt82_71)))

# Settling on glm4
N <- with(nhefs.final,length(seqn))



# Analysis restricted to N=1679
# with non-missing values on qsmk, wtchange, age, sex, and race

glm.final <- glm(qsmk~age+factor(sex)+factor(race),
                 data=nhefs.final,family="binomial")
summary(glm.final)

nhefs.final$propensity <- predict(glm.final, type="response")



require(ggplot2)
qplot(propensity,data=nhefs.final,group=qsmk,colour=qsmk,geom="density")

#####################################################
# Estimating IP weights
# Data from NHEFS
#####################################################

# Estimation of ip weights via a logistic model
fit <- glm(qsmk ~ as.factor(sex) + as.factor(race) + age + I(age^2) + 
  as.factor(education.code) + smokeintensity + 
  I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + 
  as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), 
           family = binomial(), data = nhefs0)
summary(fit)

p.qsmk.obs <- ifelse(nhefs0$qsmk == 0, 1 - predict(fit, type = "response"),
                     predict(fit, type = "response"))
nhefs0$w <- 1/p.qsmk.obs
summary(nhefs0$w)
sd(nhefs0$w)

# Estimates from a GEE
#install.packages("geepack")
require(geepack)
gee.obj <- geeglm(wt82_71~qsmk, data = nhefs0, std.err = 'san.se',
                  weights = w, id=seqn, corstr="independence")
summary(gee.obj)

library(multcomp)
#install.packages("BSagri")
library(BSagri)
 
glm.obj <- glm(wt82_71 ~ qsmk + cluster(seqn), data = nhefs0, weights = w)
summary(glm.obj)
comp<-glht(glm.obj)
CIadj<-CIGLM(comp,method="Adj")
CIadj
vcov(glm.obj)

# install.packages("sandwich")
require(sandwich)
beta <- coef(glm.obj)
SE <-sqrt(diag(vcovHC(glm.obj, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
round(cbind(beta, lcl, ucl),1)[2,] 
 
 
#####################################################
# Estimating stabilized IP weights
# Data from NHEFS
#####################################################

# estimation of denominator of ip weights
denom.fit <- glm(qsmk ~ as.factor(sex) + as.factor(race) + age + I(age^2) + 
  as.factor(education.code) + smokeintensity + 
  I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + 
  as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), 
                 family = binomial(), data = nhefs0)
summary(denom.fit)

denom.p <- predict(denom.fit, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(qsmk~1, family = binomial(), data = nhefs0)
summary(numer.fit)
numer.p <- predict(numer.fit, type = "response")

nhefs0$sw <- ifelse(nhefs0$qsmk == 0, ((1-numer.p)/(1-denom.p)),
                    (numer.p/denom.p))

summary(nhefs0$sw)

library(multcomp)
#install.packages("BSagri")
library(BSagri)

glm.obj <- glm(wt82_71~ qsmk + cluster(seqn), data = nhefs0, weights = sw)
summary(glm.obj)
comp<-glht(glm.obj)
CIadj<-CIGLM(comp,method="Adj")
CIadj

beta <- coef(glm.obj)
SE <-sqrt(diag(vcovHC(glm.obj, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
round(cbind(beta, lcl, ucl),1)[2,]

#####################################################
# Estimating the parameters of a marginal structural mean model
# with a continuous treatment Data from NHEFS
#####################################################

# Analysis restricted to subjects reporting <=25 cig/day at baseline
nhefs1 <- subset(nhefs0, smokeintensity <=25)
dim(nhefs1)

# estimation of denominator of ip weights
den.fit.obj <- lm(smkintensity82_71 ~ as.factor(sex) + 
  as.factor(race) + age + I(age^2) + 
  as.factor(education.code) + smokeintensity + I(smokeintensity^2) +
  smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71  + 
  I(wt71^2), data = nhefs1)
p.den <- predict(den.fit.obj, type = "response")
dens.den <- dnorm(nhefs1$smkintensity82_71, p.den, summary(den.fit.obj)$sigma)

# estimation of numerator of ip weights
num.fit.obj <- lm(smkintensity82_71 ~ 1, data = nhefs1)
p.num <- predict(num.fit.obj, type = "response")
dens.num <- dnorm(nhefs1$smkintensity82_71, p.num, summary(num.fit.obj)$sigma)

# estimation of Stabilized weights
nhefs1$sw.a = dens.num/dens.den
summary(nhefs1$sw.a)

gee.obj <- geeglm(wt82_71~smkintensity82_71 + I(smkintensity82_71^2), 
                  data = nhefs1, std.err = 'san.se', weights = sw.a, id=seqn, 
                  corstr="independence")
summary(gee.obj)
comp<-glht(gee.obj)
CIadj<-CIGLM(comp,method="Adj")
CIadj

#####################################################
# Estimating the parameters of a marginal structural logistic model
# Data from NHEFS
#####################################################

# using the dataset nhefs0 
# weights sw are also calculated 
# Estimating the parameters of a marginal structural logistic model
glm.obj <- glm(death ~ qsmk + cluster(seqn), data = nhefs0, 
               weights = sw, family = binomial())
summary(glm.obj)
comp<-glht(glm.obj)
CIadj<-CIGLM(comp,method="Adj")
UnlogCI(CIadj)
vcov(glm.obj)

beta <- coef(glm.obj)
SE <-sqrt(diag(vcovHC(glm.obj, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
round(exp(cbind(beta, lcl, ucl)),1)[2,]

#####################################################
# Assessing effect modification by sex using a marginal structural mean model
# Data from NHEFS
#####################################################

table(nhefs0$sex)

# estimation of denominator of ip weights
denom.fit <- glm(qsmk ~ as.factor(sex) + as.factor(race) + age + I(age^2) + 
  as.factor(education.code) + smokeintensity + 
  I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) + 
  as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2), 
                 family = binomial(), data = nhefs0)
summary(denom.fit)

denom.p <- predict(denom.fit, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(qsmk~as.factor(sex), family = binomial(), data = nhefs0)
summary(numer.fit)
numer.p <- predict(numer.fit, type = "response")

nhefs0$sw2 <- ifelse(nhefs0$qsmk == 0, ((1-numer.p)/(1-denom.p)),
                     (numer.p/denom.p))

summary(nhefs0$sw2)
sd(nhefs0$sw2)

# Estimating parameters of a marginal structural mean model
glm.obj <- glm(wt82_71~as.factor(qsmk) + as.factor(sex) + 
  as.factor(qsmk):as.factor(sex) + cluster(seqn), data = nhefs0, 
               weights = sw2)
summary(glm.obj)
comp<-glht(glm.obj)
CIadj<-CIGLM(comp,method="Adj")
CIadj
vcov(glm.obj)

beta <- coef(glm.obj)
SE <-sqrt(diag(vcovHC(glm.obj, type="HC0"))) # robust standard errors
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
p.value <- 2*(1-pnorm(abs(beta/SE)))
round(cbind(beta, SE, lcl, ucl, p.value),1)[5,]

