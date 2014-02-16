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



# (b) Stratification: classify individuals in 5 strata of approximately 
# equal size, what is the estimated effect of smoking cessation on 
# weight gain? How about 10 strata?


ApplyQuantiles <- function(data,num) {
  cut(data, breaks=c(quantile(data, probs = seq(0, 1, by = 1/num))), 
      labels=c(1:num),include.lowest=TRUE)
}

nhefs.final$rank5 <- with(nhefs.final,factor(ApplyQuantiles(propensity,5)))
nhefs.final$rank10 <- with(nhefs.final,factor(ApplyQuantiles(propensity,10)))


lm.rank5 <- lm(wt82_71~qsmk+rank5,data=nhefs.final)
lm.rank10 <- lm(wt82_71~qsmk+rank10,data=nhefs.final)

summary(lm.rank5)
summary(lm.rank10)