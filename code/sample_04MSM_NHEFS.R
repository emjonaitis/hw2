
#####################################################
# From PROGRAM 12.1 in H&R
# Descriptive statistics from NHEFS data (Table 12.1)
#####################################################
 
# reading external data file
#install.packages("sas7bdat")
library(sas7bdat)
nhefs <- read.sas7bdat("nhefs_book.sas7bdat")
nhefs <- as.data.frame(nhefs)

nhefs$cens <- as.numeric(is.na(nhefs$wt82))
nhefs$older <- as.numeric(nhefs$age > 50 & !is.na(nhefs$age))
#install.packages("car")
library(car)
nhefs$education.code <- recode(nhefs$school, " 0:8 = 1; 9:11 = 2; 12 = 3; 
                         13:15 = 4; 16:hi = 5; NA = 6 ")
nhefs$education <- recode(nhefs$education.code, " 1 = '1. 8th grade or less';
                    2 = '2. HS dropout';
                          3 = '3. HS';
                          4 = '4. College dropout';
                          5 = '5. College or more';
                          6 = 'Unknown' ")

# Analysis restricted to N=1566 
# with non-missing values in the following covariates
nhefs.original <- nhefs # Original data saved for later use
nhefs$id <- 1:nrow(nhefs)
nhefs2 <- nhefs[c("id", "qsmk", "sex", "race", "age", "school", 
                  "smokeintensity", "smokeyrs", "exercise", "active", 
                  "wt71", "wt82")]
dim(nhefs2)
nhefs2 <- as.data.frame(na.omit(nhefs2))
dim(nhefs2)
nhefs <- subset(nhefs, id %in% nhefs2$id)
rm(nhefs2)
dim(nhefs)

# mean weight change in those with and without smoking cessation
summary(nhefs$wt82_71[nhefs$cens == 0 & nhefs$qsmk == 0])
sd(nhefs$wt82_71[nhefs$cens == 0 & nhefs$qsmk == 0])
summary(nhefs$wt82_71[nhefs$cens == 0 & nhefs$qsmk == 1])
sd(nhefs$wt82_71[nhefs$cens == 0 & nhefs$qsmk == 1])

summary(glm(wt82_71 ~ qsmk, data = subset(nhefs, cens == 0)))

# restricting data for uncensored
nhefs0 <- subset(nhefs, cens == 0)

# baseline characteristics
years1 <-mean(nhefs0$age[nhefs0$qsmk == 1])
years0 <-mean(nhefs0$age[nhefs0$qsmk == 0])
male1 <-100*mean(nhefs0$sex[nhefs0$qsmk == 1]==0)
male0 <-100*mean(nhefs0$sex[nhefs0$qsmk == 0]==0)
white1 <-100*mean(nhefs0$race[nhefs0$qsmk == 1]==0)
white0 <-100*mean(nhefs0$race[nhefs0$qsmk == 0]==0)
university1 <-100*mean(nhefs0$education.code[nhefs0$qsmk == 1]==5)
university0 <-100*mean(nhefs0$education.code[nhefs0$qsmk == 0]==5)
kg1 <-mean(nhefs0$wt71[nhefs0$qsmk == 1])
kg0 <-mean(nhefs0$wt71[nhefs0$qsmk == 0])
cigs1 <-mean(nhefs0$smokeintensity[nhefs0$qsmk == 1])
cigs0 <-mean(nhefs0$smokeintensity[nhefs0$qsmk == 0])
smoke1 <-mean(nhefs0$smokeyrs[nhefs0$qsmk == 1])
smoke0 <-mean(nhefs0$smokeyrs[nhefs0$qsmk == 0])
noexer1 <-100*mean(nhefs0$exercise[nhefs0$qsmk == 1]==2)
noexer0 <-100*mean(nhefs0$exercise[nhefs0$qsmk == 0]==2)
inactive1 <-100*mean(nhefs0$active[nhefs0$qsmk == 1]==2)
inactive0 <-100*mean(nhefs0$active[nhefs0$qsmk == 0]==2)

baseline.char <- round(matrix(c(years1, years0, male1, male0, white1, white0, 
                                university1, university0, kg1, kg0, cigs1, cigs0, smoke1, smoke0, 
                                noexer1, noexer0, inactive1, inactive0),
                              ncol = 2, byrow=T), 1) 
dimnames(baseline.char) <- list(c("age, years", "men, %", "white, %", 
                                  "university, %", "weight, kg", 
                                  "Chigarettes/day", "year smoking", 
                                  "little/no exercise, %", 
                                  "inactive daily life, %"),
                                c("Smoking cessation (A=1)","No smoking cessation (A=0)"))
print(baseline.char)

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

