#######################################################
###  eGFR slope before/after covid-19 vs pneumonia  ###
###                                                 ###
###            SCREAM 3  - July 2023                ###
###        PhD, Viyaasan Mahalingasivam             ###
###                                                 ###
###  eGFR slope - before/after - covid / pneumonia  ###
###       R Script: AL Faucon / Arvid Sj?lander     ###
#######################################################


library(haven)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(sandwich)# need this package for the standard errors



### Load the data
load("P:/SCREAM2/SCREAM2_Research/Anne-Laure Faucon/Program/SCREAM covid/pneucov_10.RData")
load("P:/SCREAM2/SCREAM2_Research/Anne-Laure Faucon/Program/SCREAM covid/pneucov_10_alive.RData")
load("P:/SCREAM2/SCREAM2_Research/Anne-Laure Faucon/Program/SCREAM covid/pneucov_without_egfrpost.RData")
load("P:/SCREAM2/SCREAM2_Research/Anne-Laure Faucon/Program/SCREAM covid/pneucov_without_egfrpost_alive.RData")


#convert to data frame
data <- as.data.frame(pneucov_10)
data <- as.data.frame(pneucov_10_alive)
data <- as.data.frame(pneucov_without_egfrpost)
data <- as.data.frame(pneucov_without_egfrpost_alive)

data$time[data$period==0] <- -data$time[data$period==0]
data <- data %>% filter(time !=0)


#----------------------------------------------------------------------------------------------
### AGE-SEX ADJUSTED MODEL ###

#specify covariates in the model
covariates <- c("period", "disease", "female", "agegroup")
formula_rhs <- paste0("~", paste(paste("time", covariates, sep="*"), collapse="+"), "+time*disease*period+disease*period")
formula <- paste0("egfr", formula_rhs)
fit <- lm(formula=as.formula(formula), data=data)

#create covid and pneumonia data sets for predictions - AFTER
time_1 <- c(0,2)
lt <- length(time_1)
mean_covariates <- colMeans(data[!duplicated(data$lopnr), covariates], na.rm=TRUE)
mean_covariates <- matrix(mean_covariates, nrow=lt, ncol=length(covariates), byrow=TRUE)
colnames(mean_covariates) <- covariates
covid <- pneumonia <- data.frame(mean_covariates)
covid$time <- pneumonia$time <- time_1
covid$disease <- 1
pneumonia$disease <- 0
covid$period <- pneumonia$period <- 1

#these lines take the repeated measures into account
v <- vcovCL(x=fit, cluster=data$lopnr)
m1 <- t(model.matrix(object=as.formula(formula_rhs), data=covid))
m0 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia))
b <- matrix(fit$coef)
est1 <- t(b)%*%m1
est0 <- t(b)%*%m0
se1 <- sqrt(diag(t(m1)%*%v%*%m1))
se0 <- sqrt(diag(t(m0)%*%v%*%m0))
l1 <- est1-1.96*se1
u1 <- est1+1.96*se1
l0 <- est0-1.96*se0
u0 <- est0+1.96*se0
p1_1 <- matrix(c(est1, l1, u1), nrow=lt)
p0_1 <- matrix(c(est0, l0, u0), nrow=lt)
#print(est1[2]-est1[1])
#print(est0[2]-est0[1])

#plot
#matplot(time_1, cbind(p1_1, p0_1), type="l", lty=c(1,2,2,1,2,2), col=c(1,1,1,2,2,2),
#  xlab="time (years)", ylab="eGFR", xlim=c(-2,2))
#legend(x="topright", legend=c("covid", "pneumonia"), col=c(1,2), lty=1, bty="n")



#create covid and pneumonia data sets for predictions - BEFORE
time_0 <- c(-2,0)
lt <- length(time_0)
mean_covariates <- colMeans(data[!duplicated(data$lopnr), covariates], na.rm=TRUE)
mean_covariates <- matrix(mean_covariates, nrow=lt, ncol=length(covariates), byrow=TRUE)
colnames(mean_covariates) <- covariates
covid <- pneumonia <- data.frame(mean_covariates)
covid$time <- pneumonia$time <- time_0
covid$disease <- 1
pneumonia$disease <- 0
covid$period <- pneumonia$period <- 0

#these lines take the repeated measures into account
v <- vcovCL(x=fit, cluster=data$lopnr)
m1 <- t(model.matrix(object=as.formula(formula_rhs), data=covid))
m0 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia))
b <- matrix(fit$coef)
est1 <- t(b)%*%m1
est0 <- t(b)%*%m0
se1 <- sqrt(diag(t(m1)%*%v%*%m1))
se0 <- sqrt(diag(t(m0)%*%v%*%m0))
l1 <- est1-1.96*se1
u1 <- est1+1.96*se1
l0 <- est0-1.96*se0
u0 <- est0+1.96*se0
p1_0 <- matrix(c(est1, l1, u1), nrow=lt)
p0_0 <- matrix(c(est0, l0, u0), nrow=lt)
print(est1[2]-est1[1])
print(est0[2]-est0[1])


#plot
#matplot(time_1, cbind(p1_1, p0_1), type="l", lty=c(1,2,2,1,2,2), col=c(1,1,1,2,2,2),
#  xlab="time (years)", ylab="eGFR", xlim=c(-2,2), ylim=c(50,100))
#legend(x="topright", legend=c("covid", "pneumonia"), col=c(1,2), lty=1, bty="n")
#
#matlines(time_0, cbind(p1_0, p0_0), type="l", lty=c(1,2,2,1,2,2), col=c(1,1,1,2,2,2),
#  xlab="time (years)", ylab="eGFR")


# Figure eGFR slope before/after covid/pneumonia
summary_data <- as.data.frame(matrix(nrow=8, ncol=6))
colnames(summary_data) <- c("disease", "period", "time", "estimate", "ci_lo", "ci_hi")
summary_data[,1] <- rep(c(0,0,1,1),2)
summary_data[,2] <- c(rep(0,4), rep(1,4))
summary_data[,3] <- c(time_0, time_0, time_1, time_1)
summary_data[,4] <- c(p0_0[,1], p1_0[,1], p0_1[,1], p1_1[,1])
summary_data[,5] <- c(p0_0[,2], p1_0[,2], p0_1[,2], p1_1[,2])
summary_data[,6] <- c(p0_0[,3], p1_0[,3], p0_1[,3], p1_1[,3])

summary_data$period <- as.factor(summary_data$period)
summary_data$disease <- as.factor(summary_data$disease)
summary_data


# Figure eGFR slope before/after covid/pneumonia
windows(5,5,12)
ggplot(aes(x = time, y = estimate,  colour = disease), data = summary_data) +
  #geom_point(size=1) +
  geom_line(linetype = "solid") + 
  geom_vline(xintercept=0, color="black", linetype="dashed") +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi, fill=disease), alpha=0.2, colour=NA) +
  scale_color_manual(name="", labels = c("Pneumonia", "Covid-19"),
                     values = c("#481567FF", "#287D8EFF")) + 
  scale_fill_manual(name="", labels = c("Pneumonia", "Covid-19"),
                    values = c("#481567FF", "#287D8EFF")) + 
  xlab("Time, years") + 
  ylab("eGFR, mL/min per 1.73m ") +
  scale_x_continuous(breaks = seq(-2,2,1), label = seq(-2,2,1)) + 
  scale_y_continuous(breaks = seq(20, 100, 20), label = seq(20,100,20), limits = c(20,100)) +
  labs(title = "eGFR trajectory - Age-sex adjusted", color = "",
       subtitle = "Covid-19 vs pneumonia") + 
  #annotate("text", x=-2, y=60, label="Mean difference in eGFR slope:") + 
  #annotate("text", x=-2, y=58, label= paste0(round(fit$coef[20],2), " [", round(confint(fit)[20,1],2), " ; ", round(confint(fit)[20,2],2), "]" )   )+
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=NA), 
        plot.title = element_text(hjust = 0.5, face="bold", size=11),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=9))


rm(est0, est1, l0, l1, m0, m1, lt, se0, se1, u0, u1, v, time_0, time_1,
  p0_0, p0_1, p1_0, p1_1, covariates, covid, pneumonia, 
  formula, formula_rhs, b)



#----------------
### Bareplot of eGFR slopes before/after covid/pneumonia - age-sex adjusted model
# Intra-individual change in eGFR slope

#Function for computing linear combination and 95%CI
#coef = vector of coefficients
#vcov = variance-covariance matrix
#b = linear combination
	#0: if we do not want to include the coef
	#1: if we want to consider the coef
	#-1: if we want to substract the coef

lincomb <- function(coef, vcov, b){
  est <- sum(coef*b)
  se <- sqrt(t(b)%*%vcov%*%b)
  ci <- c(est-1.96*se, est+1.96*se)
  list(est=est, se=se, ci=ci)
}

summary(fit)
coef <- fit$coef
vcov <- vcov(fit)
mean_cov <- data.frame(t(mean_covariates[1,]))


	# slope before covid
b <- c(0,1,rep(0,4),0,1
	, mean_cov$female, mean_cov$agegroup
	,0,0)
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_before_est <- out$est
covid_before_cilo <- out$est - 1.96*out$se 
covid_before_cihi <- out$est + 1.96*out$se
rm(b,out)

	# slope after covid
b <- c(0,1,rep(0,4),1,1
	, mean_cov$female, mean_cov$agegroup
	,0,1)
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_after_est <- out$est
covid_after_cilo <- out$est - 1.96*out$se 
covid_after_cihi <- out$est + 1.96*out$se
rm(b,out)


	# slope before pneumonia
b <- c(0,1,rep(0,4),0,0
	, mean_cov$female, mean_cov$agegroup
	,0,0)
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_before_est <- out$est
pneumo_before_cilo <- out$est - 1.96*out$se 
pneumo_before_cihi <- out$est + 1.96*out$se
rm(b,out)


	# slope after pneumonia
b <- c(0,1,rep(0,4),1,0
	, mean_cov$female, mean_cov$agegroup
	,0,0)
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_after_est <- out$est
pneumo_after_cilo <- out$est - 1.96*out$se 
pneumo_after_cihi <- out$est + 1.96*out$se
rm(b,out)


	# Diff before vs after covid
b <- c(0,0,rep(0,4),1,0, 0,0,0,1)
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_covid_est <- out$est
diff_covid_cilo <- out$est - 1.96*out$se 
diff_covid_cihi <- out$est + 1.96*out$se
diff_covid <- paste0(round(diff_covid_est,2), " (", round(diff_covid_cilo,2), ";", round(diff_covid_cihi,2), ")" )
rm(b,out)



	# Diff before vs after pneumonia
b <- c(0,0,rep(0,4),1,0, 0,0,0,0)
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_pneumo_est <- out$est
diff_pneumo_cilo <- out$est - 1.96*out$se 
diff_pneumo_cihi <- out$est + 1.96*out$se
diff_pneumo <- paste0(round(diff_pneumo_est,2), " (", round(diff_pneumo_cilo,2), ";", round(diff_pneumo_cihi,2), ")" )
rm(b,out)



	# Diff after covid vs after pneumonia
b <- c(0,0,rep(0,4),0,1, 0,0,0,1)  
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_a_covid_pneumo_est <- out$est
diff_a_covid_pneumo_cilo <- out$ci[1]
diff_a_covid_pneumo_cihi <- out$ci[2]
diff_a_covid_pneumo <- paste0(round(diff_a_covid_pneumo_est,2), " (", round(diff_a_covid_pneumo_cilo,2), ";", round(diff_a_covid_pneumo_cihi,2), ")" )
rm(b,out)


### Barplot intra-individual change in eGFR ### FIGURE ARTICLE

table_barplot <- as.data.frame(matrix(nrow=4, ncol=5))
colnames(table_barplot) <- c("disease", "period", "slope", "ci_lo", "ci_hi")
table_barplot[,1] <- c(1,1,0,0)
table_barplot[,2] <- c(0,1,0,1)
table_barplot[,3] <- c(covid_before_est,  covid_after_est,  pneumo_before_est,  pneumo_after_est)
table_barplot[,4] <- c(covid_before_cilo, covid_after_cilo, pneumo_before_cilo, pneumo_after_cilo)
table_barplot[,5] <- c(covid_before_cihi, covid_after_cihi, pneumo_before_cihi, pneumo_after_cihi)
table_barplot$disease <- as.factor(table_barplot$disease)
table_barplot$period <- as.factor(table_barplot$period)
table_barplot$hosp <- "overall"
table_barplot$diff <- c(rep(diff_covid,2), rep(diff_pneumo,2))
#table_barplot$perc_change <- c(rep(diff_covid_perc, 2), rep(diff_pneumo_perc, 2))
table_barplot

#write.table(table_barplot)

windows(5,5,12)
ggplot(data=table_barplot, aes(x=disease, y=slope, fill=period)) +
  geom_bar(stat = "identity", position=position_dodge() ) +
  geom_errorbar(aes(x=disease, ymin=ci_lo, ymax=ci_hi),
                position=position_dodge(0.9),
                width=0.1) +
  scale_fill_manual(values = c("darkblue", "turquoise3"), name="", labels=c("Before", "After")) +
  xlab("") +
  ylab("eGFR slope,\nmL/min/1.73m? per year") + 
  #labs(title = "Intra-individual change in eGFR\nbefore vs after infection\nAge-sex adjusted") +
  labs(title = "Overall population\nAge-sex adjusted") +
  scale_y_continuous(labels=seq(-10,5,5), breaks=seq(-10,5,5), limits = c(-10, 5)) +
  scale_x_discrete(labels=c("Pneumonia", "Covid-19") ) +
  geom_hline(yintercept = 0) +
  annotate("text", x=0.9, y=-8.8, label=paste0("Diff eGFR slope:", diff_pneumo), size=2.7) +
  #annotate("text", x=0.9, y=-9.4, label=paste0("% Change:", diff_pneumo_perc, "%"), size=2.7) +
  annotate("text", x=2.1, y=-8.8, label=paste0("Diff eGFR slope:", diff_covid), size=2.7) + 
  #annotate("text", x=2.1, y=-9.3, label=paste0("% Change:", diff_covid_perc, "%"), size=2.7) + 
  annotate("text", x=0.5, y=5, label="A", size=5, fontface = "bold") + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y = element_line(colour = "black"),
        legend.position = "bottom",
        plot.title = element_text(size=10, face = "bold", hjust=0.5) )


rm(covid_before_est, covid_before_cilo,  covid_before_cihi,
  covid_after_est,   covid_after_cilo,   covid_after_cihi,
  pneumo_before_est, pneumo_before_cilo, pneumo_before_cihi,
  pneumo_after_est,  pneumo_after_cilo,  pneumo_after_cihi,
  diff_covid_est,    diff_covid_cilo,    diff_covid_cihi,
  diff_pneumo_est,   diff_pneumo_cilo,   diff_pneumo_cihi,
  diff_a_covid_pneumo_est, diff_a_covid_pneumo_cilo, diff_a_covid_pneumo_cihi,
  lincomb, coef, fit, vcov, diff_covid, diff_covid_perc, diff_pneumo_perc, diff_pneumo, diff_a_covid_pneumo,
  mean_covariates, mean_cov)


### Check results slope
	# Figure with slopes

covid_before_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==0 & summary_data$time==-2]
covid_before_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==0 & summary_data$time==0]
(covid_before_0 - covid_before_2)/2

pneumo_before_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==0 & summary_data$time==-2]
pneumo_before_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==0 & summary_data$time==0]
(pneumo_before_0 - pneumo_before_2)/2


covid_after_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$time==2]
covid_after_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$time==0]
(covid_after_2 - covid_after_0)/2

pneumo_after_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$time==2]
pneumo_after_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$time==0]
(pneumo_after_2 - pneumo_after_0)/2

	# Barplot with difference in eGFR slope
table_barplot

rm(covid_before_0, covid_before_2, pneumo_before_0, pneumo_before_2,
   covid_after_0, covid_after_2, pneumo_after_0, pneumo_after_2)



#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------


### FULLY ADJUSTED MODEL ###


#specify covariates in the model
covariates <- c("period", "disease",
		    "female", "agegroup",  
                "income", "educ",
                "diab", "htn", "cvd", "aki", "adm_5y", "hxpneum",  "hosp", 
                "cancer",
                "slope_before", "egfr0", "pre_creat",
                "rasi", "immunosup")

formula_rhs <- paste0("~", paste(paste("time", covariates, sep="*"), collapse="+"), "+time*disease*period+disease*period")
formula <- paste0("egfr", formula_rhs)
fit <- lm(formula=as.formula(formula), data=data)
#summary(fit)
#View(fit$coef)


#create covid and pneumonia data sets for predictions - AFTER
time_1 <- c(0,2)
lt <- length(time_1)
mean_covariates <- colMeans(data[!duplicated(data$lopnr), covariates], na.rm=TRUE)
mean_covariates <- matrix(mean_covariates, nrow=lt, ncol=length(covariates), byrow=TRUE)
colnames(mean_covariates) <- covariates
covid <- pneumonia <- data.frame(mean_covariates)
covid$time <- pneumonia$time <- time_1
covid$disease <- 1
pneumonia$disease <- 0
covid$period <- pneumonia$period <- 1

#these lines take the repeated measures into account
v <- vcovCL(x=fit, cluster=data$lopnr)
m1 <- t(model.matrix(object=as.formula(formula_rhs), data=covid))
m0 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia))
b <- matrix(fit$coef)
est1 <- t(b)%*%m1
est0 <- t(b)%*%m0
se1 <- sqrt(diag(t(m1)%*%v%*%m1))
se0 <- sqrt(diag(t(m0)%*%v%*%m0))
l1 <- est1-1.96*se1
u1 <- est1+1.96*se1
l0 <- est0-1.96*se0
u0 <- est0+1.96*se0
p1_1 <- matrix(c(est1, l1, u1), nrow=lt)
p0_1 <- matrix(c(est0, l0, u0), nrow=lt)
#print(est1[2]-est1[1])
#print(est0[2]-est0[1])

#plot
#matplot(time_1, cbind(p1_1, p0_1), type="l", lty=c(1,2,2,1,2,2), col=c(1,1,1,2,2,2),
#  xlab="time (years)", ylab="eGFR", xlim=c(-2,2))
#legend(x="topright", legend=c("covid", "pneumonia"), col=c(1,2), lty=1, bty="n")



#create covid and pneumonia data sets for predictions - BEFORE
time_0 <- c(-2,0)
lt <- length(time_0)
mean_covariates <- colMeans(data[!duplicated(data$lopnr), covariates], na.rm=TRUE)
mean_covariates <- matrix(mean_covariates, nrow=lt, ncol=length(covariates), byrow=TRUE)
colnames(mean_covariates) <- covariates
covid <- pneumonia <- data.frame(mean_covariates)
covid$time <- pneumonia$time <- time_0
covid$disease <- 1
pneumonia$disease <- 0
covid$period <- pneumonia$period <- 0

#these lines take the repeated measures into account
v <- vcovCL(x=fit, cluster=data$lopnr)
m1 <- t(model.matrix(object=as.formula(formula_rhs), data=covid))
m0 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia))
b <- matrix(fit$coef)
est1 <- t(b)%*%m1
est0 <- t(b)%*%m0
se1 <- sqrt(diag(t(m1)%*%v%*%m1))
se0 <- sqrt(diag(t(m0)%*%v%*%m0))
l1 <- est1-1.96*se1
u1 <- est1+1.96*se1
l0 <- est0-1.96*se0
u0 <- est0+1.96*se0
p1_0 <- matrix(c(est1, l1, u1), nrow=lt)
p0_0 <- matrix(c(est0, l0, u0), nrow=lt)
print(est1[2]-est1[1])
print(est0[2]-est0[1])



# Figure eGFR slope before/after covid/pneumonia
summary_data <- as.data.frame(matrix(nrow=8, ncol=6))
colnames(summary_data) <- c("disease", "period", "time", "estimate", "ci_lo", "ci_hi")
summary_data[,1] <- rep(c(0,0,1,1),2)
summary_data[,2] <- c(rep(0,4), rep(1,4))
summary_data[,3] <- c(time_0, time_0, time_1, time_1)
summary_data[,4] <- c(p0_0[,1], p1_0[,1], p0_1[,1], p1_1[,1])
summary_data[,5] <- c(p0_0[,2], p1_0[,2], p0_1[,2], p1_1[,2])
summary_data[,6] <- c(p0_0[,3], p1_0[,3], p0_1[,3], p1_1[,3])

summary_data$period <- as.factor(summary_data$period)
summary_data$disease <- as.factor(summary_data$disease)
summary_data

# Figure eGFR slope before/after covid/pneumonia
windows(5,5,12)
ggplot(aes(x = time, y = estimate,  colour = disease), data = summary_data) +
  #geom_point(size=1) +
  geom_line(linetype = "solid") + 
  geom_vline(xintercept=0, color="black", linetype="dashed") +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi, fill=disease), alpha=0.2, colour=NA) +
  scale_color_manual(name="", labels = c("Pneumonia", "Covid-19"),
                     values = c("#481567FF", "#287D8EFF")) + 
  scale_fill_manual(name="", labels = c("Pneumonia", "Covid-19"),
                    values = c("#481567FF", "#287D8EFF")) + 
  xlab("Time, years") + 
  ylab("eGFR, mL/min per 1.73m ") +
  scale_x_continuous(breaks = seq(-2,2,1), label = seq(-2,2,1)) + 
  scale_y_continuous(breaks = seq(75, 95, 5), label = seq(75,95,5), limits = c(75,95)) +
  labs(title = "eGFR trajectory - Adjusted", color = "",
       subtitle = "Covid-19 vs pneumonia") + 
  #annotate("text", x=-2, y=60, label="Mean difference in eGFR slope:") + 
  #annotate("text", x=-2, y=58, label= paste0(round(fit$coef[20],2), " [", round(confint(fit)[20,1],2), " ; ", round(confint(fit)[20,2],2), "]" )   )+
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=NA), 
        plot.title = element_text(hjust = 0.5, face="bold", size=11),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=9))



# Figure eGFR slope after covid/pneumonia - Inter-individual change in eGFR slope #### FIGURE ARTICLE
windows(5,5,12)
plot_slope_all <- ggplot(aes(x = time, y = estimate,  colour = disease), data = summary_data[summary_data$period==1,]) +
  #geom_point(size=1) +
  geom_line(linetype = "solid") + 
  #geom_vline(xintercept=0, color="black", linetype="dashed") +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi, fill=disease), alpha=0.2, colour=NA) +
  scale_color_manual(name="", labels = c("Pneumonia", "Covid-19"),
                     values = c("#481567FF", "#287D8EFF")) + 
  scale_fill_manual(name="", labels = c("Pneumonia", "Covid-19"),
                    values = c("#481567FF", "#287D8EFF")) + 
  xlab("Time, years") + 
  ylab("eGFR, mL/min per 1.73m?") +
  scale_x_continuous(breaks = seq(0,2,0.5), label = seq(0,2,0.5)) +
	## !! Change the scale according to the dataset "pneucov_10" / "pneucov_10_alive" 
##  scale_y_continuous(breaks = seq(82, 86, 1), label = seq(82,86,1), limits = c(80.7,86.5)) +
  scale_y_continuous(breaks = seq(76, 88, 2), label = seq(76,88,2), limits = c(76,88.4)) + # Fig 2
#  scale_y_continuous(breaks = seq(82, 94, 2), label = seq(82,94,2), limits = c(82, 94)) +	  # fig S7
  labs(title = "eGFR slope after infection", color = "",
       subtitle = "Overall population") + 
	## !!! Run the code below to get the values
	# Fig 2
  annotate("text", x=0.72, y=77.0, label=paste0("eGFR slope after Covid: ", slope_after_covid_perc, " % per year"), size=2.7) + 		# in ml/min: covid_after
  annotate("text", x=0.72, y=76.5, label=paste0("eGFR slope after pneumonia: ", slope_after_pneumo_perc, " % per year"), size=2.7) +	# in ml/min: pneumo_after
  annotate("text", x=0.78, y=76.0, label=paste0("Delta eGFR slope covid vs pneumonia: ", diff_a_covid_pneumo, " mL/min/1.73m? per year"), size=2.7) +
	# Fig S7
#  annotate("text", x=0.72, y=83.0, label=paste0("eGFR slope after covid: ", slope_after_covid_perc, " % per year"), size=2.7) +		# covid_after
#  annotate("text", x=0.72, y=82.5, label=paste0("eGFR slope after pneumonia: ", slope_after_pneumo_perc, " % per year"), size=2.7) +	# in ml/min: pneumo_after
#  annotate("text", x=0.90, y=82.0, label=paste0("Delta eGFR slope covid vs pneumonia: ", diff_a_covid_pneumo, " mL/min/1.73m? per year"), size=2.7) +

  annotate("text", x=-0.1, y=88, label= "A", size=5, fontface=2) +	# Fig 2
#  annotate("text", x=-0.1, y=94, label= "A", size=5, fontface=2) +	# Fig S7

  theme_classic() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=NA), 
        plot.title = element_text(hjust = 0.5, face="bold", size=11),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=9))





rm(est0, est1, l0, l1, m0, m1, lt, se0, se1, u0, u1, v, time_0, time_1,
  p0_0, p0_1, p1_0, p1_1, covariates, covid, pneumonia, 
  formula, formula_rhs, b)



#----------------
### Bareplot of eGFR slopes before/after covid/pneumonia - Adjusted model
 

#Function for computing linear combination and 95%CI
#coef = vector of coefficients
#vcov = variance-covariance matrix
#b = linear combination
	#0: if we do not want to include the coef
	#1: if we want to consider the coef
	#-1: if we want to substract the coef

lincomb <- function(coef, vcov, b){
  est <- sum(coef*b)
  se <- sqrt(t(b)%*%vcov%*%b)
  ci <- c(est-1.96*se, est+1.96*se)
  list(est=est, se=se, ci=ci)
}


summary(fit)
#View(fit$coef)
coef <- fit$coef
vcov <- vcov(fit)
mean_cov <- data.frame(t(mean_covariates[1,]))


  # Adjusted slope before covid
b <- c(0,1,rep(0,20),1
       , mean_cov$female, mean_cov$agegroup, mean_cov$income, mean_cov$educ, mean_cov$diab, mean_cov$htn
       , mean_cov$cvd, mean_cov$aki, mean_cov$adm_5y, mean_cov$hxpneum, mean_cov$hosp, mean_cov$cancer, mean_cov$slope_before 
       , mean_cov$egfr0, mean_cov$pre_creat, mean_cov$rasi, mean_cov$immunosup
       ,0,0) # b1+b22 + mean(covariate)*coefficient in the right order +++
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_before_est <- out$est
covid_before_cilo <- out$ci[1]
covid_before_cihi <- out$ci[2]
rm(b,out)


  # Adjusted slope after covid
b <- c(0,1,rep(0,19),1,1
       , mean_cov$female, mean_cov$agegroup, mean_cov$income, mean_cov$educ, mean_cov$diab, mean_cov$htn
       , mean_cov$cvd, mean_cov$aki, mean_cov$adm_5y, mean_cov$hxpneum, mean_cov$hosp, mean_cov$cancer, mean_cov$slope_before 
       , mean_cov$egfr0, mean_cov$pre_creat, mean_cov$rasi, mean_cov$immunosup
       ,0,1) # b1+b21+b22+b41
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_after_est <- out$est
covid_after_cilo <- out$ci[1]
covid_after_cihi <- out$ci[2]
covid_after <- paste0(round(covid_after_est,2), " (", round(covid_after_cilo,2), ";", round(covid_after_cihi,2), ")" )
rm(b,out)


  # slope before pneumonia
b <- c(0,1,rep(0,19),0,0
       , mean_cov$female, mean_cov$agegroup, mean_cov$income, mean_cov$educ, mean_cov$diab, mean_cov$htn
       , mean_cov$cvd, mean_cov$aki, mean_cov$adm_5y, mean_cov$hxpneum, mean_cov$hosp, mean_cov$cancer, mean_cov$slope_before 
       , mean_cov$egfr0, mean_cov$pre_creat, mean_cov$rasi, mean_cov$immunosup
       ,0,0)	#b1
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_before_est <- out$est
pneumo_before_cilo <- out$ci[1]
pneumo_before_cihi <- out$ci[2]
rm(b,out)


  # slope after pneumonia
b <- c(0,1,rep(0,19),1,0
       , mean_cov$female, mean_cov$agegroup, mean_cov$income, mean_cov$educ, mean_cov$diab, mean_cov$htn
       , mean_cov$cvd, mean_cov$aki, mean_cov$adm_5y, mean_cov$hxpneum, mean_cov$hosp, mean_cov$cancer, mean_cov$slope_before 
       , mean_cov$egfr0, mean_cov$pre_creat, mean_cov$rasi, mean_cov$immunosup
       ,0,0)  # b1+b21
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_after_est <- out$est
pneumo_after_cilo <- out$ci[1]
pneumo_after_cihi <- out$ci[2]
pneumo_after <- paste0(round(pneumo_after_est,2), " (", round(pneumo_after_cilo,2), ";", round(pneumo_after_cihi,2), ")" )
rm(b,out)


 
  # Diff before vs after covid
b <- c(rep(0,21),1,0, rep(0,17),0,1) # b21+b41
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_covid_est <- out$est
diff_covid_cilo <- out$ci[1]
diff_covid_cihi <- out$ci[2]
diff_covid <- paste0(round(diff_covid_est,2), " (", round(diff_covid_cilo,2), ";", round(diff_covid_cihi,2), ")" )
rm(b,out)



  # Diff before vs after pneumonia
b <- c(rep(0,21),1,0, rep(0,17),0,0)  #b21
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_pneumo_est <- out$est
diff_pneumo_cilo <- out$ci[1]
diff_pneumo_cihi <- out$ci[2]
diff_pneumo <- paste0(round(diff_pneumo_est,2), " (", round(diff_pneumo_cilo,2), ";", round(diff_pneumo_cihi,2), ")" )
rm(b,out)



  # Diff after covid vs after pneumonia
b <- c(rep(0,21),0,1, rep(0,17),0,1)  #b22 + b41
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_a_covid_pneumo_est <- out$est
diff_a_covid_pneumo_cilo <- out$ci[1]
diff_a_covid_pneumo_cihi <- out$ci[2]
diff_a_covid_pneumo <- paste0(round(diff_a_covid_pneumo_est,2), " (", round(diff_a_covid_pneumo_cilo,2), " ; ", round(diff_a_covid_pneumo_cihi,2), ")" )
rm(b,out)



### Calculate adjusted eGFR slope after covid (in %) ; after pneumonia
# = adjusted slope after disease / adjusted intercept of the disease of interest * 100
# 95%CI are calculated using the deltamethod() function of the msm package
# 95%CI could also be calculated using bootstrap

library(msm)

# Delta Method - Theory
#fit <- lm(formula=as.formula(formula), data=data)
#b <- fit$coef
#v <- vcov(fit)
#use delta method for the desired ratio
#Only use a linear combination of Xi in the deltamethod() function
#se <- deltamethod(g=~(x1+x2)/x2, mean=b, cov=v) 


fit
b <- fit$coef
v <- vcov(fit)

# eGFR slope after covid / adjusted intercept at T0 * 100
	# eGFR slope after covid = x2+x22+x23+x42 + sum of each "mean" covariates x their respective "covariate:time" coefficient)
	# Adjusted intercept = sum("mean" of each covariates x their respective "covariate:time" coefficient)


	# Create objects with the mean values (mean age, etc.) of each covariates ; 
	# "mean_cov" cannot be read by the deltamethod() function
female 	<- mean_cov$female
agegroup 	<- mean_cov$agegroup
income 	<- mean_cov$income
educ 		<- mean_cov$educ
diab 		<- mean_cov$diab
htn 		<- mean_cov$htn
cvd 		<- mean_cov$cvd
aki 		<- mean_cov$aki
adm_5y 	<- mean_cov$adm_5y
hxpneum 	<- mean_cov$hxpneum
hosp 		<- mean_cov$hosp
cancer 	<- mean_cov$cancer
slope_before <- mean_cov$slope_before 
egfr0 	<- mean_cov$egfr0
pre_creat 	<- mean_cov$pre_creat
rasi 		<- mean_cov$rasi
immunosup 	<- mean_cov$immunosup


se <- deltamethod(g = ~ (x2 + x22 + x23 + x42 
				+ female*x24 + agegroup*x25 + income*x26 + educ*x27 + diab*x28 + htn*x29 
				+ cvd*x30 + aki*x31 + adm_5y*x32 + hxpneum*x33 + hosp*x34 + cancer*x35 
				+ slope_before*x36 + egfr0*x37 + pre_creat*x38 + rasi*x39 + immunosup*x40) /
				(x1 + x3 + x4 + female*x5 + agegroup*x6 + income*x7 + educ*x8 + diab*x9 + htn*x10 
				+ cvd*11 + aki*12 + adm_5y*13 + hxpneum*x14 + hosp*x15 + cancer*x16 
				+ slope_before*x17 + egfr0*x18 + pre_creat*x19 + rasi*x20 + immunosup*x21 + x41)
				, mean=b, cov=v)

intercept_after_covid <- summary_data[which(summary_data$disease==1 & summary_data$period==1 & summary_data$time==0), 4]

slope_after_covid_perc <- paste0(round(abs(covid_after_est) / intercept_after_covid * 100,1), " (",
					round((abs(covid_after_est) / intercept_after_covid - 1.96*se) * 100,1), " ; ",
					round((abs(covid_after_est) / intercept_after_covid + 1.96*se) * 100,1), ")" )
slope_after_covid_perc



# eGFR slope after pneumonia / adjusted intercept at T0 * 100

se <- deltamethod(g = ~ (x2 + x22 
				+ female*x24 + agegroup*x25 + income*x26 + educ*x27 + diab*x28 + htn*x29 
				+ cvd*x30 + aki*x31 + adm_5y*x32 + hxpneum*x33 + hosp*x34 + cancer*x35 
				+ slope_before*x36 + egfr0*x37 + pre_creat*x38 + rasi*x39 + immunosup*x40) /
				(x1 + x3     + female*x5 + agegroup*x6 + income*x7 + educ*x8 + diab*x9 + htn*x10 
				+ cvd*11 + aki*12 + adm_5y*13 + hxpneum*x14 + hosp*x15 + cancer*x16 
				+ slope_before*x17 + egfr0*x18 + pre_creat*x19 + rasi*x20 + immunosup*x21  )
				, mean=b, cov=v)

intercept_after_pneumo <- summary_data[which(summary_data$disease==0 & summary_data$period==1 & summary_data$time==0), 4]

slope_after_pneumo_perc <- paste0(round(abs(pneumo_after_est) / intercept_after_pneumo * 100,1), " (",
					round((abs(pneumo_after_est) / intercept_after_pneumo - 1.96*se) * 100,1), " ; ",
					round((abs(pneumo_after_est) / intercept_after_pneumo + 1.96*se) * 100,1), ")" )

slope_after_pneumo_perc 








### Barplot

table_barplot <- as.data.frame(matrix(nrow=4, ncol=5))
colnames(table_barplot) <- c("disease", "period", "slope", "ci_lo", "ci_hi")
table_barplot[,1] <- c(1,1,0,0)
table_barplot[,2] <- c(0,1,0,1)
table_barplot[,3] <- c(covid_before_est,  covid_after_est,  pneumo_before_est,  pneumo_after_est)
table_barplot[,4] <- c(covid_before_cilo, covid_after_cilo, pneumo_before_cilo, pneumo_after_cilo)
table_barplot[,5] <- c(covid_before_cihi, covid_after_cihi, pneumo_before_cihi, pneumo_after_cihi)
table_barplot$disease <- as.factor(table_barplot$disease)
table_barplot$period <- as.factor(table_barplot$period)
table_barplot$diff <- diff_a_covid_pneumo
table_barplot


summary_data_overall <- summary_data[summary_data$period==1,]
summary_data_overall$hosp <- "overall"
summary_data_overall$slope <- c(rep(pneumo_after,2), rep(covid_after,2))
summary_data_overall$diff <- diff_a_covid_pneumo
summary_data_overall

#write.table(summary_data_overall)


windows(5,5,12)
ggplot(data=table_barplot, aes(x=disease, y=slope, fill=period)) +
  geom_bar(stat = "identity", position=position_dodge() ) +
  geom_errorbar(aes(x=disease, ymin=ci_lo, ymax=ci_hi),
                position=position_dodge(0.9),
                width=0.1) +
  scale_fill_manual(values = c("darkblue", "turquoise3"), name="", labels=c("Before", "After")) +
  xlab("") +
  ylab("eGFR slope,\nmL/min/1.73m  per year") + 
  labs(title = "Change in eGFR - Adjusted\nCovid-19 vs pneumonia") +
  scale_y_continuous(labels=seq(-4,3,1), breaks=seq(-4,3,1), limits = c(-4, 3)) +
  scale_x_discrete(labels=c("Pneumonia", "Covid-19") ) +
  geom_hline(yintercept = 0) +
  annotate("text", x=1.5, y=-3.3, label=paste0("Diff eGFR slope covid: ", diff_covid, " mL/min/1.73m?/year"), size=2.7) +
  annotate("text", x=1.5, y=-3.6, label=paste0("Diff eGFR slope pneumonia: ", diff_pneumo, " mL/min/1.73m?/year"), size=2.7) +
  annotate("text", x=1.5, y=-3.9, label=paste0("Diff eGFR slope after covid vs pneumonia: ", diff_a_covid_pneumo, " mL/min/1.73m?/year"), size=2.7) +
#  annotate("text", x=1, y=-2.6, label="Diff eGFR slope pneumonia:", size=2.7) +
#  annotate("text", x=1, y=-2.9, label=diff_pneumo, size=2.7) +  
#  annotate("text", x=2, y=-2.6, label="Diff eGFR slope Covid:", size=2.7) + 
#  annotate("text", x=2, y=-2.9, label=diff_covid, size=2.7) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y = element_line(colour = "black"),
        legend.position = "bottom",
        plot.title = element_text(size=10, face = "bold", hjust=0.5) )


rm(covid_before_est, covid_before_cilo,  covid_before_cihi,
  covid_after_est,   covid_after_cilo,   covid_after_cihi,
  pneumo_before_est, pneumo_before_cilo, pneumo_before_cihi,
  pneumo_after_est,  pneumo_after_cilo,  pneumo_after_cihi,
  diff_covid_est,    diff_covid_cilo,    diff_covid_cihi,
  diff_pneumo_est,   diff_pneumo_cilo,   diff_pneumo_cihi,
  diff_a_covid_pneumo_est, diff_a_covid_pneumo_cilo, diff_a_covid_pneumo_cihi,
  lincomb, coef, fit, vcov, diff_covid, diff_pneumo, diff_a_covid_pneumo,
  mean_covariates, mean_cov)


### Check results slope -> the results are the same
	# Figure with the slopes
covid_before_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==0 & summary_data$time==-2]
covid_before_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==0 & summary_data$time==0]
(covid_before_0 - covid_before_2)/2

pneumo_before_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==0 & summary_data$time==-2]
pneumo_before_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==0 & summary_data$time==0]
(pneumo_before_0 - pneumo_before_2)/2


covid_after_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$time==2]
covid_after_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$time==0]
(covid_after_2 - covid_after_0)/2

pneumo_after_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$time==2]
pneumo_after_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$time==0]
(pneumo_after_2 - pneumo_after_0)/2

	# Barplot with difference in eGFR slope
table_barplot


rm(covid_before_0, covid_before_2, pneumo_before_0, pneumo_before_2,
   covid_after_0, covid_after_2, pneumo_after_0, pneumo_after_2)

