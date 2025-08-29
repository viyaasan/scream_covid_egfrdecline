#######################################################
###  eGFR slope before/after covid-19 vs pneumonia  ###
###                                                 ###
###            SCREAM 3  - July 2023                ###
###        PhD, Viyaasan Mahalingasivam             ###
###                                                 ###
###  eGFR slope - before/after - covid / pneumonia  ###
###		  stratified on hospitalization	    ###
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


# Create variable diseas_hosp
data$disease_hosp <- ifelse(data$disease==0 & data$hosp==0, 0,
                     ifelse(data$disease==0 & data$hosp==1, 1,
                     ifelse(data$disease==1 & data$hosp==0,2,3))) 

data$disease_hosp <- data$disease_hosp



#----------------------------------------------------------------------------------------------
### AGE-SEX ADJUSTED MODEL ###


#specify covariates in the model
covariates <- c("period", "female","agegroup")
formula_rhs <- paste0("~", paste(paste("time", c("factor(disease_hosp, levels=0:3)", covariates), sep="*"), collapse="+"), "+time*factor(disease_hosp, levels=0:3)*period+factor(disease_hosp, levels=0:3)*period")
formula <- paste0("egfr", formula_rhs)
fit <- lm(formula=as.formula(formula), data=data)

#create covid and pneumonia data sets for predictions - AFTER
time_1 <- c(0,2)
lt <- length(time_1)
mean_covariates <- colMeans(data[!duplicated(data$lopnr), covariates], na.rm=TRUE)
mean_covariates <- matrix(mean_covariates, nrow=lt, ncol=length(covariates), byrow=TRUE)
colnames(mean_covariates) <- covariates
covid_hosp <- covid_nonhosp <- pneumonia_hosp <- pneumonia_nonhosp <- data.frame(mean_covariates)
covid_hosp$time <- covid_nonhosp$time <- pneumonia_hosp$time <- pneumonia_nonhosp$time <- time_1
pneumonia_nonhosp$disease_hosp <- 0
pneumonia_hosp$disease_hosp <- 1
covid_nonhosp$disease_hosp <- 2
covid_hosp$disease_hosp <- 3
covid_hosp$period <- covid_nonhosp$period <- pneumonia_hosp$period <- pneumonia_nonhosp$period <- 1


# make predictions, these lines take the repeated measures into account
v <- vcovCL(x=fit, cluster=data$lopnr)
m0 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia_nonhosp))
m1 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia_hosp))
m2 <- t(model.matrix(object=as.formula(formula_rhs), data=covid_nonhosp))
m3 <- t(model.matrix(object=as.formula(formula_rhs), data=covid_hosp))

b <- matrix(fit$coef)
est0 <- t(b)%*%m0
est1 <- t(b)%*%m1
est2 <- t(b)%*%m2
est3 <- t(b)%*%m3
se0 <- sqrt(diag(t(m0)%*%v%*%m0))
se1 <- sqrt(diag(t(m1)%*%v%*%m1))
se2 <- sqrt(diag(t(m2)%*%v%*%m2))
se3 <- sqrt(diag(t(m3)%*%v%*%m3))

l0 <- est0-1.96*se0
u0 <- est0+1.96*se0
l1 <- est1-1.96*se1
u1 <- est1+1.96*se1
l2 <- est2-1.96*se2
u2 <- est2+1.96*se2
l3 <- est3-1.96*se3
u3 <- est3+1.96*se3

p0_1 <- matrix(c(est0, l0, u0), nrow=lt)
p1_1 <- matrix(c(est1, l1, u1), nrow=lt)
p2_1 <- matrix(c(est2, l2, u2), nrow=lt)
p3_1 <- matrix(c(est3, l3, u3), nrow=lt)




#create covid and pneumonia data sets for predictions - BEFORE
time_0 <- c(-2,0)
lt <- length(time_0)
mean_covariates <- colMeans(data[!duplicated(data$lopnr), covariates], na.rm=TRUE)
mean_covariates <- matrix(mean_covariates, nrow=lt, ncol=length(covariates), byrow=TRUE)
colnames(mean_covariates) <- covariates
covid_hosp <- covid_nonhosp <- pneumonia_hosp <- pneumonia_nonhosp <- data.frame(mean_covariates)
covid_hosp$time <- covid_nonhosp$time <- pneumonia_hosp$time <- pneumonia_nonhosp$time <- time_0
pneumonia_nonhosp$disease_hosp <- 0
pneumonia_hosp$disease_hosp <- 1
covid_nonhosp$disease_hosp <- 2
covid_hosp$disease_hosp <- 3
covid_hosp$period <- covid_nonhosp$period <- pneumonia_hosp$period <- pneumonia_nonhosp$period <- 0


# make predictions, these lines take the repeated measures into account
v <- vcovCL(x=fit, cluster=data$lopnr)
m0 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia_nonhosp))
m1 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia_hosp))
m2 <- t(model.matrix(object=as.formula(formula_rhs), data=covid_nonhosp))
m3 <- t(model.matrix(object=as.formula(formula_rhs), data=covid_hosp))

b <- matrix(fit$coef)
est0 <- t(b)%*%m0
est1 <- t(b)%*%m1
est2 <- t(b)%*%m2
est3 <- t(b)%*%m3
se0 <- sqrt(diag(t(m0)%*%v%*%m0))
se1 <- sqrt(diag(t(m1)%*%v%*%m1))
se2 <- sqrt(diag(t(m2)%*%v%*%m2))
se3 <- sqrt(diag(t(m3)%*%v%*%m3))

l0 <- est0-1.96*se0
u0 <- est0+1.96*se0
l1 <- est1-1.96*se1
u1 <- est1+1.96*se1
l2 <- est2-1.96*se2
u2 <- est2+1.96*se2
l3 <- est3-1.96*se3
u3 <- est3+1.96*se3

p0_0 <- matrix(c(est0, l0, u0), nrow=lt)
p1_0 <- matrix(c(est1, l1, u1), nrow=lt)
p2_0 <- matrix(c(est2, l2, u2), nrow=lt)
p3_0 <- matrix(c(est3, l3, u3), nrow=lt)




# Figure eGFR slope before/after covid/pneumonia, stratified by hospitalization
summary_data <- as.data.frame(matrix(nrow=16, ncol=7))
colnames(summary_data) <- c("disease", "period", "time", "hosp", "estimate", "ci_lo", "ci_hi")
summary_data[,1] <- rep(c(0,0,0,0,1,1,1,1),2)
summary_data[,2] <- c(rep(0,8), rep(1,8))
summary_data[,3] <- c(time_0, time_0, time_0, time_0, time_1, time_1, time_1, time_1)
summary_data[,4] <- rep(c(0,0,1,1,2,2,3,3),2)
summary_data[,5] <- c(p0_0[,1], p1_0[,1], p2_0[,1], p3_0[,1],   p0_1[,1], p1_1[,1], p2_1[,1], p3_1[,1])
summary_data[,6] <- c(p0_0[,2], p1_0[,2], p2_0[,2], p3_0[,2],   p0_1[,2], p1_1[,2], p2_1[,2], p3_1[,2])
summary_data[,7] <- c(p0_0[,3], p1_0[,3], p2_0[,3], p3_0[,3],   p0_1[,3], p1_1[,3], p2_1[,3], p3_1[,3])

summary_data$period <- as.factor(summary_data$period)
summary_data$disease <- as.factor(summary_data$disease)
summary_data$hosp <- as.factor(summary_data$hosp)
summary_data$hosp2 <- as.factor(ifelse(summary_data$hosp==0 | summary_data$hosp==2, 0, 1))

summary_data



# Figure eGFR slope before/after covid/pneumonia
windows(8,5,12)
ggplot(aes(x = time, y = estimate,  colour = hosp), data = summary_data) +
  geom_line(linetype = "solid") + 
  geom_vline(xintercept=0, color="black", linetype="dashed") +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi, fill=hosp), alpha=0.2, colour=NA) +
  scale_color_manual(name="", labels = c("Pneumonia, non-hosp", "Pneumonia, hosp", "Covid-19, non-hosp", "Covid-19, hosp"),
                     values = c("darkblue", "#287D8EFF", "gold", "orange")) + 
  scale_fill_manual(name="", labels = c("Pneumonia, non-hosp", "Pneumonia, hosp", "Covid-19, non-hosp", "Covid-19, hosp"),
                    values = c("darkblue", "#287D8EFF", "gold", "orange")) + 
  xlab("Time, years") + 
  ylab("eGFR, mL/min per 1.73m?") +
  scale_x_continuous(breaks = seq(-2,2,1), label = seq(-2,2,1)) + 
  scale_y_continuous(breaks = seq(65, 92, 5), label = seq(65,92,5), limits = c(65,92)) +
  labs(title = "eGFR trajectory - Age-sex adjusted", color = "",
       subtitle = "Pneumonia vs Covid, by hospitalization") + 
  #annotate("text", x=-2, y=60, label="Mean difference in eGFR slope:") + 
  #annotate("text", x=-2, y=58, label= paste0(round(fit$coef[20],2), " [", round(confint(fit)[20,1],2), " ; ", round(confint(fit)[20,2],2), "]" )   )+
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=NA), 
        plot.title = element_text(hjust = 0.5, face="bold", size=11),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=9))


rm(est0, est1, est2, est3, l0, l1, l2, l3, m0, m1, m2, m3, lt, se0, se1, se2, se3, u0, u1, u2, u3, v, time_0, time_1,
  p0_0, p0_1, p1_0, p1_1, p2_0, p2_1, p3_0, p3_1, covariates, 
  covid_hosp, pneumonia_hosp,   covid_nonhosp, pneumonia_nonhosp, 
  formula, formula_rhs, b)




#--------------------------------------------------------------------------------------
### Bareplot of eGFR slopes before/after covid/pneumonia, stratified on hospitalization


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


	# Adjusted slope before pneumonia - non hosp
b <- c(0,1,rep(0,6),0,0,0,0
      , mean_cov$female, mean_cov$agegroup
	,0,0,0,0,0,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_nonhosp_before_est <- out$est
pneumo_nonhosp_before_cilo <- out$ci[1]
pneumo_nonhosp_before_cihi <- out$ci[2]
pneumo_nonhosp_before<- paste0(round(pneumo_nonhosp_before_est,2), " (", round(pneumo_nonhosp_before_cilo,2), " ; ", round(pneumo_nonhosp_before_cihi,2), ")" )
rm(b,out)


	# Adjusted slope before pneumonia - hosp
b <- c(0,1,rep(0,6),1,0,0,0
      , mean_cov$female, mean_cov$agegroup
	,0,0,0,0,0,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_hosp_before_est <- out$est
pneumo_hosp_before_cilo <- out$ci[1]
pneumo_hosp_before_cihi <- out$ci[2]
pneumo_hosp_before<- paste0(round(pneumo_hosp_before_est,2), " (", round(pneumo_hosp_before_cilo,2), " ; ", round(pneumo_hosp_before_cihi,2), ")" )
rm(b,out)


	# Adjusted slope before covid - non hosp
b <- c(0,1,rep(0,6),0,1,0,0
      , mean_cov$female, mean_cov$agegroup
	,0,0,0,0,0,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_nonhosp_before_est <- out$est
covid_nonhosp_before_cilo <- out$ci[1]
covid_nonhosp_before_cihi <- out$ci[2]
covid_nonhosp_before<- paste0(round(covid_nonhosp_before_est,2), " (", round(covid_nonhosp_before_cilo,2), " ; ", round(covid_nonhosp_before_cihi,2), ")" )
rm(b,out)


	# Adjusted slope before covid - hosp
b <- c(0,1,rep(0,6),0,0,1,0
      , mean_cov$female, mean_cov$agegroup
	,0,0,0,0,0,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_hosp_before_est <- out$est
covid_hosp_before_cilo <- out$ci[1]
covid_hosp_before_cihi <- out$ci[2]
covid_hosp_before<- paste0(round(covid_hosp_before_est,2), " (", round(covid_hosp_before_cilo,2), " ; ", round(covid_hosp_before_cihi,2), ")" )
rm(b,out)


	# Adjusted slope after pneumonia - non hosp
b <- c(0,1,rep(0,6),0,0,0,1
       , mean_cov$female, mean_cov$agegroup
	,0,0,0,0,0,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_nonhosp_after_est <- out$est
pneumo_nonhosp_after_cilo <- out$ci[1]
pneumo_nonhosp_after_cihi <- out$ci[2]
pneumo_nonhosp_after <- paste0(round(pneumo_nonhosp_after_est,2), " (", round(pneumo_nonhosp_after_cilo,2), " ; ", round(pneumo_nonhosp_after_cihi,2), ")" )
rm(b,out)


	# Adjusted slope after pneumonia - hosp
b <- c(0,1,rep(0,6),1,0,0,1
       , mean_cov$female, mean_cov$agegroup
	,0,0,0,1,0,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_hosp_after_est <- out$est
pneumo_hosp_after_cilo <- out$ci[1]
pneumo_hosp_after_cihi <- out$ci[2]
pneumo_hosp_after <- paste0(round(pneumo_hosp_after_est,2), " (", round(pneumo_hosp_after_cilo,2), " ; ", round(pneumo_hosp_after_cihi,2), ")" )
rm(b,out)


	# Adjusted slope after covid - non hosp
b <- c(0,1,rep(0,6),0,1,0,1
       , mean_cov$female, mean_cov$agegroup
	,0,0,0,0,1,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_nonhosp_after_est <- out$est
covid_nonhosp_after_cilo <- out$ci[1]
covid_nonhosp_after_cihi <- out$ci[2]
covid_nonhosp_after <- paste0(round(covid_nonhosp_after_est,2), " (", round(covid_nonhosp_after_cilo,2), " ; ", round(covid_nonhosp_after_cihi,2), ")" )
rm(b,out)


	# Adjusted slope after covid - hosp
b <- c(0,1,rep(0,6),0,0,1,1
       , mean_cov$female, mean_cov$agegroup
	,0,0,0,0,0,1) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_hosp_after_est <- out$est
covid_hosp_after_cilo <- out$ci[1]
covid_hosp_after_cihi <- out$ci[2]
covid_hosp_after <- paste0(round(covid_hosp_after_est,2), " (", round(covid_hosp_after_cilo,2), " ; ", round(covid_hosp_after_cihi,2), ")" )
rm(b,out)


	# Diff before vs after pneumonia - non hosp
b <- c(0,0,rep(0,6),0,0,0,1, 0,0,  0,0,0,0,0,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_pneumo_nonhosp_est <- out$est
diff_pneumo_nonhosp_cilo <- out$ci[1]
diff_pneumo_nonhosp_cihi <- out$ci[2]
diff_pneumo_nonhosp <- paste0(round(diff_pneumo_nonhosp_est,2), " (", round(diff_pneumo_nonhosp_cilo,2), " ; ", round(diff_pneumo_nonhosp_cihi,2), ")" )
rm(b,out)



	# Diff before vs after pneumonia - hosp
b <- c(0,0,rep(0,6),0,0,0,1, 0,0  ,0,0,0,1,0,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_pneumo_hosp_est <- out$est
diff_pneumo_hosp_cilo <- out$ci[1]
diff_pneumo_hosp_cihi <- out$ci[2]
diff_pneumo_hosp <- paste0(round(diff_pneumo_hosp_est,2), " (", round(diff_pneumo_hosp_cilo,2), " ; ", round(diff_pneumo_hosp_cihi,2), ")" )
rm(b,out)


	# Diff before vs after covid	- non hosp
b <- c(0,0,rep(0,6),0,0,0,1, 0,0 ,0,0,0,0,1,0)
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_covid_nonhosp_est <- out$est
diff_covid_nonhosp_cilo <- out$ci[1]
diff_covid_nonhosp_cihi <- out$ci[2]
diff_covid_nonhosp <- paste0(round(diff_covid_nonhosp_est,2), " (", round(diff_covid_nonhosp_cilo,2), " ; ", round(diff_covid_nonhosp_cihi,2), ")" )
rm(b,out)



	# Diff before vs after covid	- hosp
b <- c(0,0,rep(0,6),0,0,0,1, 0,0  ,0,0,0,0,0,1) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_covid_hosp_est <- out$est
diff_covid_hosp_cilo <- out$ci[1]
diff_covid_hosp_cihi <- out$ci[2]
diff_covid_hosp <- paste0(round(diff_covid_hosp_est,2), " (", round(diff_covid_hosp_cilo,2), " ; ", round(diff_covid_hosp_cihi,2), ")" )
rm(b,out)





### Barplot intra-individual change in eGFR ### FIGURE ARTICLE

table_barplot <- as.data.frame(matrix(nrow=8, ncol=6))
colnames(table_barplot) <- c("disease", "hosp", "period", "slope", "ci_lo", "ci_hi")
table_barplot[,1] <- c(1,1,0,0, 1,1,0,0)
table_barplot[,2] <- c(0,1,0,1, 0,1,0,1)
table_barplot[,3] <- c(0,0,0,0, 1,1,1,1)

table_barplot[,4] <- c(covid_nonhosp_before_est,  covid_hosp_before_est,  pneumo_nonhosp_before_est,  pneumo_hosp_before_est,
				covid_nonhosp_after_est,  covid_hosp_after_est,  pneumo_nonhosp_after_est,  pneumo_hosp_after_est)
table_barplot[,5] <- c(covid_nonhosp_before_cilo, covid_hosp_before_cilo, pneumo_nonhosp_before_cilo, pneumo_hosp_before_cilo,
				covid_nonhosp_after_cilo, covid_hosp_after_cilo, pneumo_nonhosp_after_cilo, pneumo_hosp_after_cilo)
table_barplot[,6] <- c(covid_nonhosp_before_cihi, covid_hosp_before_cihi, pneumo_nonhosp_before_cihi, pneumo_hosp_before_cihi,
				covid_nonhosp_after_cihi, covid_hosp_after_cihi, pneumo_nonhosp_after_cihi, pneumo_hosp_after_cihi)
table_barplot$disease<- as.factor(table_barplot$disease)
table_barplot$hosp   <- as.factor(table_barplot$hosp)
table_barplot$period <- as.factor(table_barplot$period)
table_barplot$diff   <- rep(c(diff_covid_nonhosp, diff_covid_hosp, diff_pneumo_nonhosp, diff_pneumo_hosp), 2)
#table_barplot$perc_change <- rep(c(diff_covid_nonhosp_perc, diff_covid_hosp_perc, diff_pneumo_nonhosp_perc, diff_pneumo_hosp_perc), 2)
table_barplot

#write.table(table_barplot)


	# NON HOSPITALIZED
#windows(5,5,12)
plot_nonhosp <- ggplot(data=table_barplot[table_barplot$hosp==0, ], aes(x=disease, y=slope, fill=period)) +
  geom_bar(stat = "identity", position=position_dodge() ) +
  geom_errorbar(aes(x=disease, ymin=ci_lo, ymax=ci_hi),
                position=position_dodge(0.9),
                width=0.1) +
#  scale_fill_manual(values = c("orange", "red3"), name="", labels=c("Before", "After")) +
  scale_fill_manual(values = c("darkblue", "turquoise3"), name="", labels=c("Before", "After")) +
  xlab("") +
  ylab("eGFR slope,\nmL/min/1.73m? per year") + 
#  labs(title = "Adjusted eGFR slopes - Age-sex adjusted\nBefore vs after infection\nNon hospitalized") +
  labs(title = "Non-hospitalized\nAge-sex adjusted") +
  scale_y_continuous(labels=seq(-10,5,5), breaks=seq(-10,5,5), limits = c(-10, 5)) +
  scale_x_discrete(labels=c("Pneumonia", "Covid-19") ) +
  geom_hline(yintercept = 0) +
#  annotate("text", x=1.0, y=3.6, label=paste0("Diff eGFR slope: ", diff_pneumo_nonhosp), size=2.6) +
#  annotate("text", x=1.0, y=3.0, label=paste0("% Change: ", diff_pneumo_nonhosp_perc, " %"), size=2.6) +
#  annotate("text", x=2.0, y=3.6, label=paste0("Diff eGFR slope: ", diff_covid_nonhosp), size=2.6) +
#  annotate("text", x=2.0, y=3.0, label=paste0("% Change: ", diff_covid_nonhosp_perc, " %"), size=2.6) +
  annotate("text", x=0.9, y=-8.8, label=paste0("Diff eGFR slope: ", diff_pneumo_nonhosp), size=2.6) +
  #annotate("text", x=0.9, y=-9.4, label=paste0("% Change: ", diff_pneumo_nonhosp_perc, " %"), size=2.6) +
  annotate("text", x=2.1, y=-8.8, label=paste0("Diff eGFR slope: ", diff_covid_nonhosp), size=2.6) +
  #annotate("text", x=2.1, y=-9.4, label=paste0("% Change: ", diff_covid_nonhosp_perc, " %"), size=2.6) +
  annotate("text", x=0.5, y=5, label="B", size=5, fontface = "bold") + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y = element_line(colour = "black"),
        legend.position = "bottom",
        plot.title = element_text(size=10, face = "bold", hjust=0.5) )



	# HOSPITALIZED
windows(5,5,12)
plot_hosp <- ggplot(data=table_barplot[table_barplot$hosp==1, ], aes(x=disease, y=slope, fill=period)) +
  geom_bar(stat = "identity", position=position_dodge() ) +
  geom_errorbar(aes(x=disease, ymin=ci_lo, ymax=ci_hi),
                position=position_dodge(0.9),
                width=0.1) +
#  scale_fill_manual(values = c("orange", "red3"), name="", labels=c("Before", "After")) +
  scale_fill_manual(values = c("darkblue", "turquoise3"), name="", labels=c("Before", "After")) +
  xlab("") +
  ylab("eGFR slope,\nmL/min/1.73m? per year") + 
#  labs(title = "Adjusted eGFR slopes - Age-sex adjusted\nBefore vs after infection\nHospitalized") +
  labs(title = "Hospitalized\nAge-sex adjusted") +
  scale_y_continuous(labels=seq(-10,5,5), breaks=seq(-10,5,5), limits = c(-10, 5)) +
  scale_x_discrete(labels=c("Pneumonia", "Covid-19") ) +
  geom_hline(yintercept = 0) +
#  annotate("text", x=1.0, y=3.6, label=paste0("Diff eGFR slope: ", diff_pneumo_hosp), size=2.6) +
#  annotate("text", x=1.0, y=3.0, label=paste0("% Change: ", diff_pneumo_hosp_perc, " %"), size=2.6) +
#  annotate("text", x=2.0, y=3.6, label=paste0("Diff eGFR slope: ", diff_covid_hosp), size=2.6) +
#  annotate("text", x=2.0, y=3.0, label=paste0("% Change: ", diff_covid_hosp_perc, " %"), size=2.6) +
  annotate("text", x=0.9, y=-8.8, label=paste0("Diff eGFR slope: ", diff_pneumo_hosp), size=2.6) +
  #annotate("text", x=0.9, y=-9.4, label=paste0("% Change: ", diff_pneumo_hosp_perc, " %"), size=2.6) +
  annotate("text", x=2.1, y=-8.8, label=paste0("Diff eGFR slope: ", diff_covid_hosp), size=2.6) +
  #annotate("text", x=2.1, y=-9.4, label=paste0("% Change: ", diff_covid_hosp_perc, " %"), size=2.6) +
  annotate("text", x=0.5, y=5, label="C", size=5, fontface = "bold") + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y = element_line(colour = "black"),
        legend.position = "bottom",
        plot.title = element_text(size=10, face = "bold", hjust=0.5) )



windows(9,5,12)
grid.arrange(plot_nonhosp, plot_hosp, ncol = 2)


rm(covid_nonhosp_before_est, 	covid_nonhosp_before_cilo,	covid_nonhosp_before_cihi, 	covid_nonhosp_before,
  covid_hosp_before_est,     	covid_hosp_before_cilo,     	covid_hosp_before_cihi, 	covid_hosp_before,
  pneumo_nonhosp_before_est, 	pneumo_nonhosp_before_cilo, 	pneumo_nonhosp_before_cihi, 	pneumo_nonhosp_before,
  pneumo_hosp_before_est,    	pneumo_hosp_before_cilo,    	pneumo_hosp_before_cihi, 	pneumo_hosp_before,
  covid_nonhosp_after_est,   	covid_nonhosp_after_cilo,   	covid_nonhosp_after_cihi, 	covid_nonhosp_after,
  covid_hosp_after_est,      	covid_hosp_after_cilo,      	covid_hosp_after_cihi, 		covid_hosp_after,
  pneumo_nonhosp_after_est,  	pneumo_nonhosp_after_cilo,  	pneumo_nonhosp_after_cihi, 	pneumo_nonhosp_after,
  pneumo_hosp_after_est,     	pneumo_hosp_after_cilo,     	pneumo_hosp_after_cihi, 	pneumo_hosp_after,
  diff_covid_hosp_est,		diff_covid_hosp_cilo, 		diff_covid_hosp_cihi,		diff_covid_hosp, 	diff_covid_hosp_perc,
  diff_covid_nonhosp_est,	diff_covid_nonhosp_cilo, 	diff_covid_nonhosp_cihi,	diff_covid_nonhosp, diff_covid_nonhosp_perc,
  diff_pneumo_hosp_est,		diff_pneumo_hosp_cilo, 		diff_pneumo_hosp_cihi,		diff_pneumo_hosp,  diff_pneumo_hosp_perc,
  diff_pneumo_nonhosp_est,	diff_pneumo_nonhosp_cilo, 	diff_pneumo_nonhosp_cihi,	diff_pneumo_nonhosp, diff_pneumo_nonhosp_perc,
  lincomb, coef, fit, vcov, mean_covariates, mean_cov,
  plot_hosp, plot_nonhosp)



### Check results slope -> the results are the same
	# Figure with the slopes - before
pneumo_nonhosp_before_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==0 & summary_data$hosp==0 & summary_data$time==-2]
pneumo_nonhosp_before_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==0 & summary_data$hosp==0 & summary_data$time==0]
(pneumo_nonhosp_before_0 - pneumo_nonhosp_before_2)/2

pneumo_hosp_before_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==0 & summary_data$hosp==1 & summary_data$time==-2]
pneumo_hosp_before_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==0 & summary_data$hosp==1 & summary_data$time==0]
(pneumo_hosp_before_0 - pneumo_hosp_before_2)/2

covid_nonhosp_before_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==0 & summary_data$hosp==2 & summary_data$time==-2]
covid_nonhosp_before_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==0 & summary_data$hosp==2 & summary_data$time==0]
(covid_nonhosp_before_0 - covid_nonhosp_before_2)/2

covid_hosp_before_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==0 & summary_data$hosp==3 & summary_data$time==-2]
covid_hosp_before_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==0 & summary_data$hosp==3 & summary_data$time==0]
(covid_hosp_before_0 - covid_hosp_before_2)/2

	# Figure with the slopes - after
pneumo_nonhosp_after_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$hosp==0 & summary_data$time==2]
pneumo_nonhosp_after_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$hosp==0 & summary_data$time==0]
(pneumo_nonhosp_after_2 - pneumo_nonhosp_after_0)/2

pneumo_hosp_after_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$hosp==1 & summary_data$time==2]
pneumo_hosp_after_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$hosp==1 & summary_data$time==0]
(pneumo_hosp_after_2 - pneumo_hosp_after_0)/2

covid_nonhosp_after_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$hosp==2 & summary_data$time==2]
covid_nonhosp_after_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$hosp==2 & summary_data$time==0]
(covid_nonhosp_after_2 - covid_nonhosp_after_0)/2

covid_hosp_after_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$hosp==3 & summary_data$time==2]
covid_hosp_after_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$hosp==3 & summary_data$time==0]
(covid_hosp_after_2 - covid_hosp_after_0)/2

	# Barplot with difference in eGFR slope
table_barplot


rm(covid_nonhosp_before_0, covid_nonhosp_before_2, pneumo_nonhosp_before_0, pneumo_nonhosp_before_2,
   covid_hosp_before_0, covid_hosp_before_2, pneumo_hosp_before_0, pneumo_hosp_before_2,
   covid_nonhosp_after_0, covid_nonhosp_after_2, pneumo_nonhosp_after_0, pneumo_nonhosp_after_2,
   covid_hosp_after_0, covid_hosp_after_2, pneumo_hosp_after_0, pneumo_hosp_after_2)




#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
### ADJUSTED MODEL ###

#specify covariates in the model
covariates <- c("period", "female","agegroup",  
                "income", "educ",
                "diab", "htn", "cvd", "aki", "adm_5y", "hxpneum",  
                "cancer",
                "slope_before", "egfr0", "pre_creat",
                "rasi", "immunosup")
formula_rhs <- paste0("~", paste(paste("time", c("factor(disease_hosp, levels=0:3)", covariates), sep="*"), collapse="+"), "+time*factor(disease_hosp, levels=0:3)*period+factor(disease_hosp, levels=0:3)*period")
formula <- paste0("egfr", formula_rhs)
fit <- lm(formula=as.formula(formula), data=data)

#create covid and pneumonia data sets for predictions - AFTER
time_1 <- c(0,2)
lt <- length(time_1)
mean_covariates <- colMeans(data[!duplicated(data$lopnr), covariates], na.rm=TRUE)
mean_covariates <- matrix(mean_covariates, nrow=lt, ncol=length(covariates), byrow=TRUE)
colnames(mean_covariates) <- covariates
covid_hosp <- covid_nonhosp <- pneumonia_hosp <- pneumonia_nonhosp <- data.frame(mean_covariates)
covid_hosp$time <- covid_nonhosp$time <- pneumonia_hosp$time <- pneumonia_nonhosp$time <- time_1
pneumonia_nonhosp$disease_hosp <- 0
pneumonia_hosp$disease_hosp <- 1
covid_nonhosp$disease_hosp <- 2
covid_hosp$disease_hosp <- 3
covid_hosp$period <- covid_nonhosp$period <- pneumonia_hosp$period <- pneumonia_nonhosp$period <- 1


# make predictions, these lines take the repeated measures into account
v <- vcovCL(x=fit, cluster=data$lopnr)
m0 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia_nonhosp))
m1 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia_hosp))
m2 <- t(model.matrix(object=as.formula(formula_rhs), data=covid_nonhosp))
m3 <- t(model.matrix(object=as.formula(formula_rhs), data=covid_hosp))

b <- matrix(fit$coef)
est0 <- t(b)%*%m0
est1 <- t(b)%*%m1
est2 <- t(b)%*%m2
est3 <- t(b)%*%m3
se0 <- sqrt(diag(t(m0)%*%v%*%m0))
se1 <- sqrt(diag(t(m1)%*%v%*%m1))
se2 <- sqrt(diag(t(m2)%*%v%*%m2))
se3 <- sqrt(diag(t(m3)%*%v%*%m3))

l0 <- est0-1.96*se0
u0 <- est0+1.96*se0
l1 <- est1-1.96*se1
u1 <- est1+1.96*se1
l2 <- est2-1.96*se2
u2 <- est2+1.96*se2
l3 <- est3-1.96*se3
u3 <- est3+1.96*se3

p0_1 <- matrix(c(est0, l0, u0), nrow=lt)
p1_1 <- matrix(c(est1, l1, u1), nrow=lt)
p2_1 <- matrix(c(est2, l2, u2), nrow=lt)
p3_1 <- matrix(c(est3, l3, u3), nrow=lt)



#create covid and pneumonia data sets for predictions - BEFORE
time_0 <- c(-2,0)
lt <- length(time_0)
mean_covariates <- colMeans(data[!duplicated(data$lopnr), covariates], na.rm=TRUE)
mean_covariates <- matrix(mean_covariates, nrow=lt, ncol=length(covariates), byrow=TRUE)
colnames(mean_covariates) <- covariates
covid_hosp <- covid_nonhosp <- pneumonia_hosp <- pneumonia_nonhosp <- data.frame(mean_covariates)
covid_hosp$time <- covid_nonhosp$time <- pneumonia_hosp$time <- pneumonia_nonhosp$time <- time_0
pneumonia_nonhosp$disease_hosp <- 0
pneumonia_hosp$disease_hosp <- 1
covid_nonhosp$disease_hosp <- 2
covid_hosp$disease_hosp <- 3
covid_hosp$period <- covid_nonhosp$period <- pneumonia_hosp$period <- pneumonia_nonhosp$period <- 0


# make predictions, these lines take the repeated measures into account
v <- vcovCL(x=fit, cluster=data$lopnr)
m0 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia_nonhosp))
m1 <- t(model.matrix(object=as.formula(formula_rhs), data=pneumonia_hosp))
m2 <- t(model.matrix(object=as.formula(formula_rhs), data=covid_nonhosp))
m3 <- t(model.matrix(object=as.formula(formula_rhs), data=covid_hosp))

b <- matrix(fit$coef)
est0 <- t(b)%*%m0
est1 <- t(b)%*%m1
est2 <- t(b)%*%m2
est3 <- t(b)%*%m3
se0 <- sqrt(diag(t(m0)%*%v%*%m0))
se1 <- sqrt(diag(t(m1)%*%v%*%m1))
se2 <- sqrt(diag(t(m2)%*%v%*%m2))
se3 <- sqrt(diag(t(m3)%*%v%*%m3))

l0 <- est0-1.96*se0
u0 <- est0+1.96*se0
l1 <- est1-1.96*se1
u1 <- est1+1.96*se1
l2 <- est2-1.96*se2
u2 <- est2+1.96*se2
l3 <- est3-1.96*se3
u3 <- est3+1.96*se3

p0_0 <- matrix(c(est0, l0, u0), nrow=lt)
p1_0 <- matrix(c(est1, l1, u1), nrow=lt)
p2_0 <- matrix(c(est2, l2, u2), nrow=lt)
p3_0 <- matrix(c(est3, l3, u3), nrow=lt)




# Figure eGFR slope before/after covid/pneumonia, stratified by hospitalization
summary_data <- as.data.frame(matrix(nrow=16, ncol=7))
colnames(summary_data) <- c("disease", "period", "time", "hosp", "estimate", "ci_lo", "ci_hi")
summary_data[,1] <- rep(c(0,0,0,0,1,1,1,1),2)
summary_data[,2] <- c(rep(0,8), rep(1,8))
summary_data[,3] <- c(time_0, time_0, time_0, time_0, time_1, time_1, time_1, time_1)
summary_data[,4] <- rep(c(0,0,1,1,2,2,3,3),2)
summary_data[,5] <- c(p0_0[,1], p1_0[,1], p2_0[,1], p3_0[,1],   p0_1[,1], p1_1[,1], p2_1[,1], p3_1[,1])
summary_data[,6] <- c(p0_0[,2], p1_0[,2], p2_0[,2], p3_0[,2],   p0_1[,2], p1_1[,2], p2_1[,2], p3_1[,2])
summary_data[,7] <- c(p0_0[,3], p1_0[,3], p2_0[,3], p3_0[,3],   p0_1[,3], p1_1[,3], p2_1[,3], p3_1[,3])

summary_data$period <- as.factor(summary_data$period)
summary_data$disease <- as.factor(summary_data$disease)
summary_data$hosp <- as.factor(summary_data$hosp)
summary_data$hosp2 <- as.factor(ifelse(summary_data$hosp==0 | summary_data$hosp==2, 0, 1))

summary_data


# Figure eGFR slope before/after covid/pneumonia
windows(7,5,12)
ggplot(aes(x = time, y = estimate,  colour = hosp), data = summary_data) +
  geom_line(linetype = "solid") + 
  geom_vline(xintercept=0, color="black", linetype="dashed") +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi, fill=hosp), alpha=0.2, colour=NA) +
  scale_color_manual(name="", labels = c("Pneumonia, non-hosp", "Pneumonia, hosp", "Covid-19, non-hosp", "Covid-19, hosp"),
                     values = c("darkblue", "#287D8EFF", "gold", "orange")) + 
  scale_fill_manual(name="", labels = c("Pneumonia, non-hosp", "Pneumonia, hosp", "Covid-19, non-hosp", "Covid-19, hosp"),
                    values = c("darkblue", "#287D8EFF", "gold", "orange")) + 
  xlab("Time, years") + 
  ylab("eGFR, mL/min per 1.73m?") +
  scale_x_continuous(breaks = seq(-2,2,1), label = seq(-2,2,1)) + 
  scale_y_continuous(breaks = seq(75, 90, 5), label = seq(75,90,5), limits = c(75,90)) +
  labs(title = "eGFR trajectory - Adjusted", color = "",
       subtitle = "Pneumonia vs Covid, by hospitalization") + 
  #annotate("text", x=-2, y=60, label="Mean difference in eGFR slope:") + 
  #annotate("text", x=-2, y=58, label= paste0(round(fit$coef[20],2), " [", round(confint(fit)[20,1],2), " ; ", round(confint(fit)[20,2],2), "]" )   )+
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=NA), 
        plot.title = element_text(hjust = 0.5, face="bold", size=11),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=9))




# Figure eGFR slope AFTER covid/pneumonia, post infection - FIGURE ARTICLE

	# NON-HOSPITALIZED

windows(5,5,12)
plot_nonhosp <- ggplot(aes(x = time, y = estimate,  colour = hosp), data = summary_data[summary_data$period==1 & summary_data$hosp2==0, ]) +
  geom_line(linetype = "solid") + 
  #geom_vline(xintercept=0, color="black", linetype="dashed") +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi, fill=hosp), alpha=0.2, colour=NA) +
  scale_color_manual(name="", labels = c("Pneumonia", "Covid-19"),
                     values = c("darkblue", "#287D8EFF")) + 
  scale_fill_manual(name="", labels = c("Pneumonia", "Covid-19"),
                    values = c("darkblue",  "#287D8EFF")) + 
  xlab("Time, years") + 
  ylab("eGFR, mL/min per 1.73m?") +
  scale_x_continuous(breaks = seq(0,2,0.5), label = seq(0,2,0.5)) + 
#  scale_y_continuous(breaks = seq(76, 88, 2), label = seq(76,88,2), limits = c(76,88.4)) +	# Fig 2
  scale_y_continuous(breaks = seq(82, 94, 2), label = seq(82,94,2), limits = c(82, 94)) +		# fig S7

  labs(title = "eGFR slope after infection", color = "",
       subtitle = "Non-hospitalized") + 
	## !!! Run the code below (diff eGFR)
	# Fig 2
#  annotate("text", x=0.80, y=77, label=paste0("eGFR slope after covid: ", slope_after_covid_nonhosp_perc, " % per year"), size=2.7 ) + 		# in mL/min: covid_nonhosp_after
#  annotate("text", x=0.80, y=76.5, label=paste0("eGFR slope after pneumonia: ", slope_after_pneumo_nonhosp_perc, " % per year"), size=2.7 ) + 	# in ml/min: pneumo_nonhosp_after
#  annotate("text", x=0.90, y=76.0, label=paste0("Delta eGFR slope covid vs pneumonia: ", diff_covid_pneumo_nonhosp, " mL/min/1.73m? per year"), size=2.7 ) + 
	# Fig S7
  annotate("text", x=0.72, y=83.0, label=paste0("eGFR slope after covid: ", slope_after_covid_nonhosp_perc, " % per year"), size=2.7) +		# covid_after
  annotate("text", x=0.72, y=82.5, label=paste0("eGFR slope after pneumonia: ", slope_after_pneumo_nonhosp_perc, " % per year"), size=2.7) +	# in ml/min: pneumo_after
  annotate("text", x=0.90, y=82.0, label=paste0("Delta eGFR slope covid vs pneumonia: ", diff_covid_pneumo_nonhosp, " mL/min/1.73m? per year"), size=2.7) +

#  annotate("text", x=-0.1, y=88, label= "B", size=5, fontface=2) +	# Fig 2
  annotate("text", x=-0.1, y=94, label= "B", size=5, fontface=2) +	# Fig S7

  theme_classic() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=NA), 
        plot.title = element_text(hjust = 0.5, face="bold", size=11),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=9))


	# HOSPITALIZED
plot_hosp <- ggplot(aes(x = time, y = estimate,  colour = hosp), data = summary_data[summary_data$period==1 & summary_data$hosp2==1, ]) +
  geom_line(linetype = "solid") + 
  #geom_vline(xintercept=0, color="black", linetype="dashed") +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi, fill=hosp), alpha=0.2, colour=NA) +
  scale_color_manual(name="", labels = c("Pneumonia", "Covid-19"),
                     values = c("darkblue", "#287D8EFF")) + 
  scale_fill_manual(name="", labels = c("Pneumonia", "Covid-19"),
                    values = c("darkblue",  "#287D8EFF")) + 
  xlab("Time, years") + 
  ylab("eGFR, mL/min per 1.73m?") +
  scale_x_continuous(breaks = seq(0,2,0.5), label = seq(0,2,0.5)) + 
#  scale_y_continuous(breaks = seq(76, 88, 2), label = seq(76,88,2), limits = c(76,88.4)) +	# Fig 2
  scale_y_continuous(breaks = seq(82, 94, 2), label = seq(82,94,2), limits = c(82, 94)) +		# fig S7

  labs(title = "eGFR slope after infection", color = "",
       subtitle = "Hospitalized") + 
	## !!! Run the code below (diff eGFR)
	# Fig 2
#  annotate("text", x=0.8, y=77, label=paste0("eGFR slope after covid: ", slope_after_covid_hosp_perc, " % per year"), size=2.7 ) + 			# in ml/min: covid_hosp_after
#  annotate("text", x=0.8, y=76.5, label=paste0("eGFR slope after pneumonia: ", slope_after_pneumo_hosp_perc, " % per year"), size=2.7 ) + 	# in ml/min: pneumo_hosp_after
#  annotate("text", x=0.9, y=76.0, label=paste0("Delta eGFR slope covid vs pneumonia: ", diff_covid_pneumo_hosp, " mL/min/1.73m? per year"), size=2.7 ) + 
	# Fig S7
  annotate("text", x=0.72, y=83.0, label=paste0("eGFR slope after covid: ", slope_after_covid_hosp_perc, " % per year"), size=2.7) +		# covid_after
  annotate("text", x=0.72, y=82.5, label=paste0("eGFR slope after pneumonia: ", slope_after_pneumo_hosp_perc, " % per year"), size=2.7) +	# in ml/min: pneumo_after
  annotate("text", x=0.90, y=82.0, label=paste0("Delta eGFR slope covid vs pneumonia: ", diff_covid_pneumo_hosp, " mL/min/1.73m? per year"), size=2.7) +

#  annotate("text", x=-0.1, y=88, label= "C", size=5, fontface=2) +	# Fig 2
  annotate("text", x=-0.1, y=94, label= "C", size=5, fontface=2) +	# Fig S7
 
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill=NA), 
        plot.title = element_text(hjust = 0.5, face="bold", size=11),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=9))



rm(est0, est1, est2, est3, l0, l1, l2, l3, m0, m1, m2, m3, lt, se0, se1, se2, se3, u0, u1, u2, u3, v, time_0, time_1,
  p0_0, p0_1, p1_0, p1_1, p2_0, p2_1, p3_0, p3_1, covariates, 
  covid_hosp, pneumonia_hosp,   covid_nonhosp, pneumonia_nonhosp, 
  formula, formula_rhs, b, plot_hosp, plot_nonhosp)




#--------------------------------------------------------------------------------------
### Bareplot of eGFR slopes before/after covid/pneumonia, stratified on hospitalization


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



	# Adjusted slope after pneumonia - non hosp
b <- c(0,1,rep(0,20),0,0,0,1
       , mean_cov$female, mean_cov$agegroup, mean_cov$income, mean_cov$educ, mean_cov$diab, mean_cov$htn
       , mean_cov$cvd, mean_cov$aki, mean_cov$adm_5y, mean_cov$hxpneum, mean_cov$cancer, mean_cov$slope_before 
       , mean_cov$egfr0, mean_cov$pre_creat, mean_cov$rasi, mean_cov$immunosup
	,0,0,0,0,0,0) # b1+b25+ all the covariates*time
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_nonhosp_after_est <- out$est
pneumo_nonhosp_after_cilo <- out$ci[1]
pneumo_nonhosp_after_cihi <- out$ci[2]
pneumo_nonhosp_after <- paste0(round(pneumo_nonhosp_after_est,2), " (", round(pneumo_nonhosp_after_cilo,2), " ; ", round(pneumo_nonhosp_after_cihi,2), ")" )
rm(b,out)


	# Adjusted slope after pneumonia - hosp
b <- c(0,1,rep(0,20),1,0,0,1
       , mean_cov$female, mean_cov$agegroup, mean_cov$income, mean_cov$educ, mean_cov$diab, mean_cov$htn
       , mean_cov$cvd, mean_cov$aki, mean_cov$adm_5y, mean_cov$hxpneum, mean_cov$cancer, mean_cov$slope_before 
       , mean_cov$egfr0, mean_cov$pre_creat, mean_cov$rasi, mean_cov$immunosup
	,0,0,0,1,0,0) # b1+b22+b25+b45+ all the covariates*time
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
pneumo_hosp_after_est <- out$est
pneumo_hosp_after_cilo <- out$ci[1]
pneumo_hosp_after_cihi <- out$ci[2]
pneumo_hosp_after <- paste0(round(pneumo_hosp_after_est,2), " (", round(pneumo_hosp_after_cilo,2), " ; ", round(pneumo_hosp_after_cihi,2), ")" )
rm(b,out)


	# Adjusted slope after covid - non hosp
b <- c(0,1,rep(0,20),0,1,0,1
       , mean_cov$female, mean_cov$agegroup, mean_cov$income, mean_cov$educ, mean_cov$diab, mean_cov$htn
       , mean_cov$cvd, mean_cov$aki, mean_cov$adm_5y, mean_cov$hxpneum, mean_cov$cancer, mean_cov$slope_before 
       , mean_cov$egfr0, mean_cov$pre_creat, mean_cov$rasi, mean_cov$immunosup
	,0,0,0,0,1,0) # b1+b23+b25+b46+ all the covariates*time
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_nonhosp_after_est <- out$est
covid_nonhosp_after_cilo <- out$ci[1]
covid_nonhosp_after_cihi <- out$ci[2]
covid_nonhosp_after <- paste0(round(covid_nonhosp_after_est,2), " (", round(covid_nonhosp_after_cilo,2), " ; ", round(covid_nonhosp_after_cihi,2), ")" )
rm(b,out)


	# Adjusted slope after covid - hosp
b <- c(0,1,rep(0,20),0,0,1,1
       , mean_cov$female, mean_cov$agegroup, mean_cov$income, mean_cov$educ, mean_cov$diab, mean_cov$htn
       , mean_cov$cvd, mean_cov$aki, mean_cov$adm_5y, mean_cov$hxpneum, mean_cov$cancer, mean_cov$slope_before 
       , mean_cov$egfr0, mean_cov$pre_creat, mean_cov$rasi, mean_cov$immunosup
	,0,0,0,0,0,1) # b1+b24+b25+b47+ all the covariates*time
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
covid_hosp_after_est <- out$est
covid_hosp_after_cilo <- out$ci[1]
covid_hosp_after_cihi <- out$ci[2]
covid_hosp_after <- paste0(round(covid_hosp_after_est,2), " (", round(covid_hosp_after_cilo,2), " ; ", round(covid_hosp_after_cihi,2), ")" )
rm(b,out)


	# Diff before vs after pneumonia - non hosp
b <- c(0,0,rep(0,20),0,0,0,1, rep(0,16),	0,0,0,0,0,0)
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_pneumo_nonhosp_est <- out$est
diff_pneumo_nonhosp_cilo <- out$ci[1]
diff_pneumo_nonhosp_cihi <- out$ci[2]
diff_pneumo_nonhosp <- paste0(round(diff_pneumo_nonhosp_est,2), " (", round(diff_pneumo_nonhosp_cilo,2), " ; ", round(diff_pneumo_nonhosp_cihi,2), ")" )
rm(b,out)



	# Diff before vs after pneumonia - hosp
b <- c(0,0,rep(0,20),0,0,0,1, rep(0,16),	0,0,0,1,0,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_pneumo_hosp_est <- out$est
diff_pneumo_hosp_cilo <- out$ci[1]
diff_pneumo_hosp_cihi <- out$ci[2]
diff_pneumo_hosp <- paste0(round(diff_pneumo_hosp_est,2), " (", round(diff_pneumo_hosp_cilo,2), " ; ", round(diff_pneumo_hosp_cihi,2), ")" )
rm(b,out)



	# Diff before vs after covid	- non hosp
b <- c(0,0,rep(0,20),0,0,0,1, rep(0,16), 0,0,0,0,1,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_covid_nonhosp_est <- out$est
diff_covid_nonhosp_cilo <- out$ci[1]
diff_covid_nonhosp_cihi <- out$ci[2]
diff_covid_nonhosp <- paste0(round(diff_covid_nonhosp_est,2), " (", round(diff_covid_nonhosp_cilo,2), " ; ", round(diff_covid_nonhosp_cihi,2), ")" )
rm(b,out)



	# Diff before vs after covid	- hosp
b <- c(0,0,rep(0,20),0,0,0,1, rep(0,16), 0,0,0,0,0,1) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_covid_hosp_est <- out$est
diff_covid_hosp_cilo <- out$ci[1]
diff_covid_hosp_cihi <- out$ci[2]
diff_covid_hosp <- paste0(round(diff_covid_hosp_est,2), " (", round(diff_covid_hosp_cilo,2), " ; ", round(diff_covid_hosp_cihi,2), ")" )
rm(b,out)



	# Diff after covid vs after pneumonia - nonhosp
b <- c(0,0,rep(0,20),0,1,0,0, rep(0,16), 0,0,0,0,1,0) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_covid_pneumo_nonhosp_est <- out$est
diff_covid_pneumo_nonhosp_cilo <- out$ci[1]
diff_covid_pneumo_nonhosp_cihi <- out$ci[2]
diff_covid_pneumo_nonhosp <- paste0(round(diff_covid_pneumo_nonhosp_est,2), " (", round(diff_covid_pneumo_nonhosp_cilo,2), " ; ", round(diff_covid_pneumo_nonhosp_cihi,2), ")" )
rm(b,out)



	# Diff after covid vs after pneumonia - hosp
	# !! we have to put "-1", otherwise the result is not correct
b <- c(0,0,rep(0,20),-1,0,1,0, rep(0,16), 0,0,0,-1,0,1) 
out <- lincomb(coef=coef, vcov=vcov, b=b)
print(out)
diff_covid_pneumo_hosp_est <- out$est
diff_covid_pneumo_hosp_cilo <- out$ci[1]
diff_covid_pneumo_hosp_cihi <- out$ci[2]
diff_covid_pneumo_hosp <- paste0(round(diff_covid_pneumo_hosp_est,2), " (", round(diff_covid_pneumo_hosp_cilo,2), " ; ", round(diff_covid_pneumo_hosp_cihi,2), ")" )
rm(b,out)



### Calculate adjusted eGFR slope after covid (in %), after pneumonia, stratified on hospitalization
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
cancer 	<- mean_cov$cancer
slope_before <- mean_cov$slope_before 
egfr0 	<- mean_cov$egfr0
pre_creat 	<- mean_cov$pre_creat
rasi 		<- mean_cov$rasi
immunosup 	<- mean_cov$immunosup



# Delta method - after covid non hosp
se <- deltamethod(g = ~ (x2 + x24 + x26 + x47 
				+ female*x27 + agegroup*x28 + income*x29 + educ*x30 + diab*x31 + htn*x32 
				+ cvd*x33 + aki*x34 + adm_5y*x35 + hxpneum*x36 + cancer*x37 
				+ slope_before*x38 + egfr0*x39 + pre_creat*x40 + rasi*x41 + immunosup*x42) /
				(x1 + x4 + x6 + female*x7 + agegroup*x8 + income*x9 + educ*x10 + diab*x11 + htn*x12 
				+ cvd*13 + aki*14 + adm_5y*15 + hxpneum*x16 + cancer*x17 
				+ slope_before*x18 + egfr0*x19 + pre_creat*x20 + rasi*x21 + immunosup*x22 + x44)
				, mean=b, cov=v)

intercept_after_covid_nonhosp <- summary_data[which(summary_data$disease==1 & summary_data$hosp==2 & summary_data$period==1 & summary_data$time==0), 5]

slope_after_covid_nonhosp_perc <- paste0(round(abs(covid_nonhosp_after_est) / intercept_after_covid_nonhosp * 100,1), " (",
					round((abs(covid_nonhosp_after_est) / intercept_after_covid_nonhosp - 1.96*se) * 100,1), " ; ",
					round((abs(covid_nonhosp_after_est) / intercept_after_covid_nonhosp + 1.96*se) * 100,1), ")" )
slope_after_covid_nonhosp_perc




# Delta method - covid hosp
se <- deltamethod(g = ~ (x2 + x25 + x26 + x48 
				+ female*x27 + agegroup*x28 + income*x29 + educ*x30 + diab*x31 + htn*x32 
				+ cvd*x33 + aki*x34 + adm_5y*x35 + hxpneum*x36 + cancer*x37 
				+ slope_before*x38 + egfr0*x39 + pre_creat*x40 + rasi*x41 + immunosup*x42) /
				(x1 + x4 + x6 + female*x7 + agegroup*x8 + income*x9 + educ*x10 + diab*x11 + htn*x12 
				+ cvd*13 + aki*14 + adm_5y*15 + hxpneum*x16 + cancer*x17 
				+ slope_before*x18 + egfr0*x19 + pre_creat*x20 + rasi*x21 + immunosup*x22 + x45)
				, mean=b, cov=v)


intercept_after_covid_hosp <- summary_data[which(summary_data$disease==1 & summary_data$hosp==3 & summary_data$period==1 & summary_data$time==0), 5]

slope_after_covid_hosp_perc <- paste0(round(abs(covid_hosp_after_est) / intercept_after_covid_hosp * 100,1), " (",
					round((abs(covid_hosp_after_est) / intercept_after_covid_hosp - 1.96*se) * 100,1), " ; ",
					round((abs(covid_hosp_after_est) / intercept_after_covid_hosp + 1.96*se) * 100,1), ")" )
slope_after_covid_hosp_perc



# Delta method - after pneumo non hosp
se <- deltamethod(g = ~ (x2 + x26 
				+ female*x27 + agegroup*x28 + income*x29 + educ*x30 + diab*x31 + htn*x32 
				+ cvd*x33 + aki*x34 + adm_5y*x35 + hxpneum*x36 + cancer*x37 
				+ slope_before*x38 + egfr0*x39 + pre_creat*x40 + rasi*x41 + immunosup*x42) /
				(x1 + x4 + x6 + female*x7 + agegroup*x8 + income*x9 + educ*x10 + diab*x11 + htn*x12 
				+ cvd*13 + aki*14 + adm_5y*15 + hxpneum*x16 + cancer*x17 
				+ slope_before*x18 + egfr0*x19 + pre_creat*x20 + rasi*x21 + immunosup*x22   )
				, mean=b, cov=v)


intercept_after_pneumo_nonhosp <- summary_data[which(summary_data$disease==0 & summary_data$hosp==0 & summary_data$period==1 & summary_data$time==0), 5]

slope_after_pneumo_nonhosp_perc <- paste0(round(abs(pneumo_nonhosp_after_est) / intercept_after_pneumo_nonhosp * 100,1), " (",
					round((abs(pneumo_nonhosp_after_est) / intercept_after_pneumo_nonhosp - 1.96*se) * 100,1), " ; ",
					round((abs(pneumo_nonhosp_after_est) / intercept_after_pneumo_nonhosp + 1.96*se) * 100,1), ")" )
slope_after_pneumo_nonhosp_perc



# Delta method - pneumonia hosp
se <- deltamethod(g = ~ (x2 + x23 + x26 + x46 
				+ female*x27 + agegroup*x28 + income*x29 + educ*x30 + diab*x31 + htn*x32 
				+ cvd*x33 + aki*x34 + adm_5y*x35 + hxpneum*x36 + cancer*x37 
				+ slope_before*x38 + egfr0*x39 + pre_creat*x40 + rasi*x41 + immunosup*x42) /
				(x1 + x4 + x6 + female*x7 + agegroup*x8 + income*x9 + educ*x10 + diab*x11 + htn*x12 
				+ cvd*13 + aki*14 + adm_5y*15 + hxpneum*x16 + cancer*x17 
				+ slope_before*x18 + egfr0*x19 + pre_creat*x20 + rasi*x21 + immunosup*x22 + x43)
				, mean=b, cov=v)


intercept_after_pneumo_hosp <- summary_data[which(summary_data$disease==0 & summary_data$hosp==1 & summary_data$period==1 & summary_data$time==0), 5]

slope_after_pneumo_hosp_perc <- paste0(round(abs(pneumo_hosp_after_est) / intercept_after_pneumo_hosp * 100,1), " (",
					round((abs(pneumo_hosp_after_est) / intercept_after_pneumo_hosp - 1.96*se) * 100,1), " ; ",
					round((abs(pneumo_hosp_after_est) / intercept_after_pneumo_hosp + 1.96*se) * 100,1), ")" )
slope_after_pneumo_hosp_perc





### Barplot

table_barplot <- as.data.frame(matrix(nrow=4, ncol=5))
colnames(table_barplot) <- c("disease", "hosp", "slope", "ci_lo", "ci_hi")
table_barplot[,1] <- c(1,1,0,0)
table_barplot[,2] <- c(0,1,0,1)
table_barplot[,3] <- c(covid_nonhosp_after_est,  covid_hosp_after_est,  pneumo_nonhosp_after_est,  pneumo_hosp_after_est)
table_barplot[,4] <- c(covid_nonhosp_after_cilo, covid_hosp_after_cilo, pneumo_nonhosp_after_cilo, pneumo_hosp_after_cilo)
table_barplot[,5] <- c(covid_nonhosp_after_cihi, covid_hosp_after_cihi, pneumo_nonhosp_after_cihi, pneumo_hosp_after_cihi)
table_barplot$disease <- as.factor(table_barplot$disease)
table_barplot$hosp <- as.factor(table_barplot$hosp)
table_barplot



summary_data_hosp <- summary_data[summary_data$period==1,]
summary_data_hosp$slope <- c(rep(pneumo_nonhosp_after,2), rep(pneumo_hosp_after,2), 
					rep(covid_nonhosp_after,2), rep(covid_hosp_after,2) )
summary_data_hosp$diff <- rep(c(rep(diff_covid_pneumo_nonhosp,2), rep(diff_covid_pneumo_hosp,2)),2)
summary_data_hosp

#write.table(summary_data_hosp)


windows(5,5,12)
ggplot(data=table_barplot, aes(x=disease, y=slope, fill=hosp)) +
  geom_bar(stat = "identity", position=position_dodge() ) +
  geom_errorbar(aes(x=disease, ymin=ci_lo, ymax=ci_hi),
                position=position_dodge(0.9),
                width=0.1) +
  scale_fill_manual(values = c("orange", "red3"), name="", labels=c("Non-hospitalized", "Hospitalized")) +
  xlab("") +
  ylab("eGFR slope,\nmL/min/1.73m? per year") + 
  labs(title = "Adjusted eGFR slopes\nafter Covid-19 vs after pneumonia\nby hospitalization status") +
  scale_y_continuous(labels=seq(-4,3,1), breaks=seq(-4,3,1), limits = c(-4, 3)) +
  scale_x_discrete(labels=c("After pneumonia", "After Covid-19") ) +
  geom_hline(yintercept = 0) +
  annotate("text", x=1.5, y=2.3, label=paste0("Pneumonia, non-hosp: ", pneumo_nonhosp_after, " mL/min/1.73m?/y"), size=2.7) +
  annotate("text", x=1.5, y=2.0, label=paste0("Pneumonia, hosp: ",     pneumo_hosp_after, " mL/min/1.73m?/y"), size=2.7) +
  annotate("text", x=1.5, y=1.7, label=paste0("Covid-19, non-hosp: ",  covid_nonhosp_after, " mL/min/1.73m?/y"), size=2.7) +
  annotate("text", x=1.5, y=1.4, label=paste0("Covid-19, hosp: ",      covid_hosp_after, " mL/min/1.73m?/y"), size=2.7) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y = element_line(colour = "black"),
        legend.position = "bottom",
        plot.title = element_text(size=10, face = "bold", hjust=0.5) )


rm(covid_nonhosp_after_est, covid_nonhosp_after_cilo,  covid_nonhosp_after_cihi, covid_nonhosp_after,
  covid_hosp_after_est,   covid_hosp_after_cilo,   covid_hosp_after_cihi, covid_hosp_after,
  pneumo_nonhosp_after_est, pneumo_nonhosp_after_cilo, pneumo_nonhosp_after_cihi, pneumo_nonhosp_after,
  pneumo_hosp_after_est,  pneumo_hosp_after_cilo,  pneumo_hosp_after_cihi, pneumo_hosp_after,
  lincomb, coef, fit, vcov, mean_covariates, mean_cov)



### Check results slope -> the results are the same
	# Figure with the slopes
pneumo_nonhosp_after_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$hosp==0 & summary_data$time==2]
pneumo_nonhosp_after_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$hosp==0 & summary_data$time==0]
(pneumo_nonhosp_after_2 - pneumo_nonhosp_after_0)/2

pneumo_hosp_after_2 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$hosp==1 & summary_data$time==2]
pneumo_hosp_after_0 <- summary_data$estimate[summary_data$disease==0 & summary_data$period==1 & summary_data$hosp==1 & summary_data$time==0]
(pneumo_hosp_after_2 - pneumo_hosp_after_0)/2

covid_nonhosp_after_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$hosp==2 & summary_data$time==2]
covid_nonhosp_after_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$hosp==2 & summary_data$time==0]
(covid_nonhosp_after_2 - covid_nonhosp_after_0)/2

covid_hosp_after_2 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$hosp==3 & summary_data$time==2]
covid_hosp_after_0 <- summary_data$estimate[summary_data$disease==1 & summary_data$period==1 & summary_data$hosp==3 & summary_data$time==0]
(covid_hosp_after_2 - covid_hosp_after_0)/2

	# Barplot with difference in eGFR slope
table_barplot


rm(covid_nonhosp_after_0, covid_nonhosp_after_2, pneumo_nonhosp_after_0, pneumo_nonhosp_after_2,
   covid_hosp_after_0, covid_hosp_after_2, pneumo_hosp_after_0, pneumo_hosp_after_2)




