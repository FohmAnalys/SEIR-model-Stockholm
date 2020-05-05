


# Remove workspace
rm(list = ls(all=TRUE))

heading <- "
#------------------------------------------------------------------------------------
#
# FILE:         Estimate_SEIR
# CONTACT:      Lisa Brouwers
# EMAIL:        analysenheten@folkhalsomyndigheten.se
# AFFILIATION:  FD-AN
# PROJECT:      Modelling COVID-19
#
#
# Created:      2020-03-15 
# Updated:      2020-05-05
# R version:    3.5.2
#
# What the script does: This script estimates the parameters of the infectivity described in the report 
#                       'Skattning av peakdagen och antal infekterade i covid-19-utbrottet i Stockholms län'.
#                       With these we are able to estimate the number of infectious individuals at different time points, 
#                       the cumulative number of infected, and the estimated effective basic reproduction number. 
#                       If you want to reproduce the results obtained, or change anything, first write the  
#                       directory where the project folder is saved on your own computer in poject.path below.
#                       
#                       
#                       
#------------------------------------------------------------------------------------
\n"





#-----------------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------------

# input your own path where folder (end with "/") is located in project.path.
# e.g. project.path 	     <- "C:/Users/Modelling/"

project.path 	     <- "G:/Projekt/Modellering nCoV/Estimation_SEIR_model/Sharing/"


data.path 		     <- paste(project.path, "Data", sep="")
script.path 	     <- paste(project.path, "Script", sep="")
table.path         <- paste(project.path, "Results/Tables", sep="")
figure.path        <- paste(project.path, "Results/Figures", sep="")


#---------------------------------------------------------------------------------------------------
# LIBRARIES
#---------------------------------------------------------------------------------------------------
library(reshape2)
library(openxlsx)     # to write tables in excel
library(RColorBrewer)
require(rootSolve)    # to load function multiroot that finds roots to system of equations
require(deSolve)


#---------------------------------------------------------------------------------------------------
# Basic functions and set.up
#---------------------------------------------------------------------------------------------------

Blues <- brewer.pal(n = 9, name = "Blues")

roundUp <- function(x, level = c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * level[[which(x <= 10^floor(log10(x)) * level)[[1]]]]
}



## From a sample x, calculate a credibility interval with wanted level. 
CRI <- function(x, level = 0.95){
  n <- length(x)
  L = (1 - level) / 2
  U = 1 - (1 - level) / 2
  x <- sort(x)
  resL <- x[n * L] 
  resU <- x[n * U]
  #names(resL) <- paste((L * 100), "%") 
  #names(resU) <- paste((U * 100), "%")
  
  return(c(resL,resU))
  
}

CRI_90_low <- function(x){
  return(CRI(x, level = 0.9)[1])
}
CRI_90_up <- function(x){
  return(CRI(x, level = 0.9)[2])
}

CRI_95_low <- function(x){
  return(CRI(x, level = 0.95)[1])
}
CRI_95_up <- function(x){
  return(CRI(x, level = 0.95)[2])
}



#---------------------------------------------------------------------------------------------------
# Read data. 
#---------------------------------------------------------------------------------------------------

# This is the data used in the analysis. It differs from some of the reported case data since we 
# removed imported cases. 

Stockholm_Data_10_april <- read.table(file=paste(data.path,"/Data_2020-04-10Ny.txt",sep=""), sep = " ", header=TRUE)


N <- 2374550



#---------------------------------------------------------------------------------------------------
# Functions for analysis
#---------------------------------------------------------------------------------------------------

## Some important notes on notations

## p_symp := fraction reported cases. Denoted p_r in report. 
## p_asymp := fraction non-reported cases. Denoted p_o in report.
## p_lower_inf := how much lower infectivity non-reported cases have. Denoted q_o in report.
## t_b := 16 march 2020
## iter := In optimisation, how many iterations/guesses that are run.
## eta := rate of leaving incubation period. Denoted rho in report.


set.seed(854789652)

eta_value    <- 1/5.1
gammaD_value <- 1/5

## Tolerance for ode and optimisation. 
## This tolerance not always needed but better to be safe even though it takes a bit longer!
## For faster analysis use different tolerance or use default values.
Atol = 1e-8
Rtol = 1e-10

Estimate_function_Stockholm_only_local <- function(p_symp = 0.5, p_lower_inf = 0.5, gammaD = gammaD_value, eta = eta_value , iter = 20){
 
  
  
  ## Daily incidence reported cases and their dates
  Incidence <- Stockholm_Data_10_april$Incidens
  Datum     <- as.Date(Stockholm_Data_10_april$Datum)
  Day       <- as.numeric(Datum) - as.numeric(as.Date("2019-12-31"))
  
  dayatmonth <- c(1,31,29,31,30,31,30,31,31,30,31,30,31)
  dayatyear <- cumsum(dayatmonth)
  Namedate <- as.Date(dayatyear, origin = "2019-12-31")
  
    
  Opt_par_names <- c("delta","epsilon","theta")
    
  ## Function to create guesses for the optimisation
  ## The range of the guesses can be changed, 
  ## these are good for the specific dates and parameter combinations of p_symp and p_lower_inf
  Guesses <- function(){ 
    
    u_d <- runif(1, 0.05, 0.6) # guess for delta 
    u_e <- runif(1,-0.6, 0)    # guess for epsilon
    u_t <- runif(1, 0, 15)     # guess for theta

    return(c(u_d, u_e, u_t))
    
  }
  
  ## The time-dependent infectivity rate 
  beta_decrease <- function(t, delta, epsilon, theta){
    
    t_b <- as.numeric(as.Date("2020-03-16")) - as.numeric(as.Date("2019-12-31")) # Mid-point for change in infectivity, the day swedes were urged to stay at home as much as possible and to work from home
    
    res <- ((1-delta)/(1+exp(epsilon*(-(t-t_b)))) + delta)* theta 
    
    return(res)
  }
  
  beta.peak.free <- beta_decrease
  
  ## The time-dependent basic reproductive number
  Basic_repr<- function(t, delta, epsilon, theta, gamma){ 
    res <- p_symp * beta.peak.free(t, delta, epsilon, theta) / gamma + (1 - p_symp) * p_lower_inf * beta.peak.free(t, delta, epsilon, theta) / gamma
    return(res)
  
  }
 
  ## The SEIR model. 
  ## Note that the rate of going from latency to infectious is denoted eta here, rho in the report
  seir.model.asymptomatics <- function(time, state, parameters) {
    # S       <- state[1] # susceptibles
    # E       <- state[2] # latent/exposed but not infectious
    # I_symp  <- state[3] # infected who get reported
    # I_asymp <- state[4] # infected who remain non-reported
    # R       <- state[5] # recovered/immune
    par <- as.list(c(state, parameters))
    with(par, {
      dS        <- -beta.peak.free(time, delta, epsilon, theta) * S * I_symp/N - p_lower_inf*beta.peak.free(time, delta, epsilon, theta) * S * I_asymp/N
      dE        <- beta.peak.free(time, delta, epsilon, theta) * S * I_symp/N + p_lower_inf*beta.peak.free(time, delta, epsilon, theta) * S * I_asymp/N - eta*E
      dI_symp   <- p_symp * eta * E      - gammaD * I_symp
      dI_asymp  <- (1 - p_symp)* eta * E - gammaD * I_asymp
      dR        <- gammaD * (I_symp + I_asymp)
      dx        <- c(dS, dE, dI_symp, dI_asymp, dR)
      list(dx)
    }
    )
  }
  
  ## assumption on initial number of infected. 
  ## In our main analysis we start with 1 infectious individual at t_0 = 17th Ferbruary. You can instead choose the commented one 
  ## if you want to try with non-reported cases as well. 
  #init <- c(S = N - Incidence[1]*(1 + (1-p_symp)/p_symp), E = 0, I_symp = Incidence[1], I_asymp = Incidence[1]*(1-p_symp)/p_symp , R = 0)
  init <- c(S = N-Incidence[1],E = 0, I_symp = Incidence[1], I_asymp = 0 , R = 0)
  
  model <- seir.model.asymptomatics
  
  RSS <- function(parameters) {
        
    names(parameters) <- Opt_par_names
    Dummy_infectivity <- beta.peak.free(t = c(0: 700), delta = parameters[1], epsilon = parameters[2],  theta = parameters[3])
    # if the infectivity is negative, throw away guess
      if(min(Dummy_infectivity) < 0 ){
      res <- 10^12
      #print("negative infectivity")
      return(res)
    }else{
      # choose tolerance atol and rtol
      #out <- ode(y = init, times = Day, func = model, parms = parameters)
      out <- ode(y = init, times = Day, func = model, parms = parameters, atol = Atol, rtol = Rtol)
      
      
          
      fit_S <- out[ , 2]
      fit_E <- out[ , 3]
      fit_I_symp <- out[ , 4]
      fit_I_asymp <- out[ , 5]
      
      fitted_incidence  <- p_symp * fit_E * eta
      
      ## For trancparacy, the old incorrect incidence was expressed as:
      #fitted_incidence  <- beta.peak.free(out[,1], delta = parameters[1], epsilon = parameters[2],  theta = parameters[3]) * fit_S * fit_I_symp/N  
      
      return(sum((Incidence - fitted_incidence)^2))
    }
  }
  


  print("Optimisation initialised")
  Guess <- Guesses()
  #conl <- list(maxit = 1000)
  conl <- list(maxit = 1000, abstol = Atol, reltol = Rtol)
  
  Opt <- 0
  Opt <- optim(Guess, RSS, control = list(conl), hessian = TRUE)

  while(Opt$convergence>0){
    
    Guess <- Guesses()
    Opt <- optim(Guess, RSS, control = list(conl), hessian = TRUE)

  }

  i = 1 
  
  while(i <= iter){
    Guess <- Guesses()
    # choose tolerance abstol and reltol
    conl <- list(maxit = 1000, parscale = c(0.000001, 1, 0.01), abstol = Atol, reltol = Rtol)
    #conl <- list(maxit = 1000, parscale = c(0.000001, 1, 0.01)) 
    #conl <- list(maxit = 1000, parscale = c(0.0001, 0.1, 0.0001)) 
    
    Opt2 <- optim(Guess, RSS, control = list(conl), hessian = TRUE)
   
    if((Opt2$convergence == 0)){
      
      if( ((i/iter)*100) %% 10 == 0){ 
        How_much_done <- (i/iter)*100
        print(paste(How_much_done,"% done", sep=""))
      }
      
      i = i + 1
      if(Opt2$value < Opt$value){Opt <- Opt2}
    }
    
  }
  
  
  Opt_par <- Opt$par
  Opt_par <- setNames(Opt$par, Opt_par_names)
  return(list(Observed_incidence = Incidence, Population_size = N, Day = Day, dayatyear = dayatyear, Namedate = Namedate, Optimisation = Opt, Infectivity = beta.peak.free, Basic_reproduction = Basic_repr, Initial_values = init, SEIR_model = model))
}




#---------------------------------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------------------------------


gammaD <- gammaD_value
eta <- eta_value


# Analysis  p_symp_use <- 0.0127
# Analysis  p_lower_inf_use <- 1, 0.55, 0.11
p_symp_use      <- 0.0127
p_asymp_use     <- 1 - p_symp_use
p_lower_inf_use <- 1


Est_par_model <- Estimate_function_Stockholm_only_local(p_symp = p_symp_use, p_lower_inf = p_lower_inf_use)


# Days of incidence
Day <- Est_par_model$Day

#N   <- Est_par_model$Population_size

dayatyear <- Est_par_model$dayatyear
Namedate  <- Est_par_model$Namedate 

#Observed incidence based on region
Observed_incidence <-  Est_par_model$Observed_incidence 
Est <- Est_par_model$Optimisation
Est
Opt_par <- setNames(Est$par, c("delta","epsilon","theta"))
Opt_par

# functions based on model scenario
Basic_repr <- Est_par_model$Basic_reproduction
beta       <- Est_par_model$Infectivity
SEIR_model <- Est_par_model$SEIR_model

# initial values based on model scenario
init       <- Est_par_model$Initial_values


RSS_value <- Est$value


H         <- Est$hessian
sigest    <- sqrt(RSS_value/(length(Observed_incidence)-3))
NeginvH2  <- solve(1/(2*sigest^2)*H)
sdParams  <- sqrt(diag(NeginvH2))
sdParams

options("scipen"=100, "digits"=4)
#default options("scipen"=0, "digits"=7)
RSS_value

#-------------------------------------------
# Look at the results, compare with prevalence 
# 27th March to 3rd April
#-------------------------------------------




t <- (Day[1]):(Day[length(Day)]+14+11) # time in days
fit <- data.frame(ode(y = init, times = t, func = SEIR_model , parms = Opt_par))
fit_S <- fit[ , 2]
fit_E <- fit[ , 3]
fit_I_symp <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_R <- fit[ , 6]
fit_I <- fit_I_symp + fit_I_asymp
fit_I_E <- fit_E + fit_I
fit_cum_inf <- N - fit_S



## The mean prevalence same days as the Hälsorapport Stockholmsstudien (27th March to 3rd April)
Infectious         <- fit_I_symp + fit_I_asymp #+ fit_E
InfectiousF <-  Infectious[40:47]
mean(InfectiousF/N)





## Look at the estimated reported cases and fitted

fitted_incidence             <- p_symp_use * fit_E * eta
# fitted_incidence_non_report  <- (1 - p_symp_use) * fit_E * eta

plot(Observed_incidence,type="o", ylab = "Reported cases", xlab ="Day")
lines(fitted_incidence,  lwd = 2, col = "blue")
legend("topleft", c("Observed reported cases", "Fitted reported cases"), lwd=c(1,2), pch = c(1,NA), lty =c(1,1), col = c("black","blue"))



## Look at the estimated infectivity and basic reproductive number
## 1 jan, 1 feb  etc
dayatyear_march_april <- c(1, 32, 61, 61 + 31, 61 + 31 + 30, 61 + 31 + 30 + 31)
NameDateMarchApril <- as.Date(dayatyear_march_april,origin ="2019-12-31")

par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 
plot(c(0:150),Basic_repr(c(0:150), delta = Est$par[1], epsilon = Est$par[2],  theta = Est$par[3]  ,gamma = gammaD),type="l", ylab="R0(t)",lwd=2, 
     main = "Estimated reproductive number",xlab="", xaxt ='n')
abline(v=Day[1], lty = 2)
abline(v=Day[length(Day)], lty = 2)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)

plot(c(0:150),beta(c(0:150), delta = Est$par[1], epsilon = Est$par[2],  theta = Est$par[3]),type="l", ylab="Infectivity",lwd=2, 
     main = "Estimated infectivity",xlab="", xaxt='n')
abline(v=Day[1], lty = 2)
abline(v=Day[length(Day)], lty = 2)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


#-------------------------------------------------
# Calculate results to save in tables
#-------------------------------------------------


## To create bootstrap CI
CI_level_05 <- 0.025

delta_high    <- qnorm(1-CI_level_05, mean = Opt_par[1], sd = sdParams[1], lower.tail = TRUE, log.p = FALSE)
epsilon_high  <- qnorm(1-CI_level_05, mean = Opt_par[2], sd = sdParams[2], lower.tail = TRUE, log.p = FALSE)
theta_high    <- qnorm(1-CI_level_05, mean = Opt_par[3], sd = sdParams[3], lower.tail = TRUE, log.p = FALSE)


delta_low     <- max(0,qnorm(CI_level_05, mean = Opt_par[1], sd = sdParams[1], lower.tail = TRUE, log.p = FALSE))
epsilon_low   <- qnorm(CI_level_05, mean = Opt_par[2], sd = sdParams[2], lower.tail = TRUE, log.p = FALSE)
theta_low     <- qnorm(CI_level_05, mean = Opt_par[3], sd = sdParams[3], lower.tail = TRUE, log.p = FALSE)

Opt_par
delta_high
delta_low

epsilon_high
epsilon_low

theta_high
theta_low


p.v        <- rnorm(1000, mean = Opt_par[1], sd = sdParams[1]) # delta.v
p.v[which(p.v<0)] <- 0
epsilon.v  <- rnorm(1000,mean = Opt_par[2], sd = sdParams[2])
theta.v    <- rnorm(1000,mean = Opt_par[3], sd = sdParams[3])


R0.v.Dag1 <- c()
R0.v.DagSista <- c()
for(i in 1:length(p.v)){
  
  R0.v.Dag1[i]      <- Basic_repr(Day[1], delta = p.v[i], epsilon = epsilon.v[i], theta = theta.v[i], gamma = gammaD) 
  R0.v.DagSista[i]  <- Basic_repr(Day[length(Day)], delta = p.v[i], epsilon = epsilon.v[i], theta = theta.v[i], gamma = gammaD) 
}



CRI(R0.v.Dag1, level= 0.95)
CRI(R0.v.DagSista, level= 0.95)

R0_low  <- c(CRI_95_low(R0.v.Dag1), CRI_95_low(R0.v.DagSista))
R0_high <- c(CRI_95_up(R0.v.Dag1), CRI_95_up(R0.v.DagSista))

R0_Mean <- Basic_repr(Day, delta = Est$par[1], epsilon = Est$par[2],  theta = Est$par[3]  ,gamma = gammaD)



R0_Mean[1]
R0_Mean[length(Day)]


#############################################
## save estimated parameters and their SE  ##
#############################################


res_param        <- c(round(c(p_asymp_use, p_lower_inf_use),digits=3), round(mean(InfectiousF / N), digits = 5), round(c(RSS_value, Est$par[1], sdParams[1], Est$par[2], sdParams[2], Est$par[3], sdParams[3]), digits = 3))
names(res_param) <- c("p_0", "q_0","27 mars - 3 april", "RSS" ,"delta", "s.e.", "epsilon", "s.e.", "theta", "s.e.")

CIp       <- paste("[",round(delta_low,digits = 3), ", ", round(delta_high, digits = 3),"]", sep="")
CIepsilon <- paste("[",round(epsilon_low,digits = 3), ", ",round(epsilon_high, digits = 3),"]", sep="")
CItheta   <- paste("[",round(theta_low,digits = 3), ", ",round(theta_high, digits = 3),"]", sep="")

CI_param <- c("", "", "", "", CIp, "", CIepsilon, "", CItheta, "")

MAT_para <- matrix(c(res_param,CI_param),ncol = 10, nrow = 2, byrow = TRUE)

df.res <- as.data.frame(MAT_para)
colnames(df.res) <- names(res_param)



XL_file_name <- paste(table.path,"/Res_para_p_non-reported_",p_asymp_use,"_infect_",p_lower_inf_use,".xlsx", sep ="")
write.xlsx(df.res, XL_file_name )




#############################################
## Save results with 31 days forecast ##
#############################################


t <- (Day[1]):(Day[length(Day)]+31) # time in days
fit         <- data.frame(ode(y = init, times = t, func = SEIR_model , parms = Opt_par))
fit_S       <- fit[ , 2]
fit_E       <- fit[ , 3]
fit_I_symp  <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_R       <- fit[ , 6]
fit_I       <- fit_I_symp + fit_I_asymp
fit_I_E     <- fit_E + fit_I
fit_cum_inf <- N - fit_S




fitted_incidence  <- p_symp_use * fit_E * eta
fitted_incidence_non_report  <- (1 - p_symp_use) * fit_E * eta


Cum_Inf_inc <- data.frame(Cumulative = fit_cum_inf, 
                          Incidence_reported = fitted_incidence, Incidence_non_reported = fitted_incidence_non_report, Datum = as.Date(t, origin = "2019-12-31"))

fit_Cum_Inf_inc <- cbind(fit, Cum_Inf_inc)


file_name <- paste(table.path,"/Raw_data_fitted_model", "_para_p_asymp", p_asymp_use, "infect", p_lower_inf_use, ".txt", sep ="")

write.table(fit_Cum_Inf_inc, file=file_name, append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE)




##################################
## Now calculate bootstrap CI's ##
##################################


fit_S.v       <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_E.v       <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_I_symp.v  <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_I_asymp.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
Fit_I.v       <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_cum_inf.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))

fitted_incidence.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
effective_reprod.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))


for(i in 1:length(p.v)){
  Opt_parDummy = c(delta = p.v[i], epsilon = epsilon.v[i], theta = theta.v[i])
  
  fitDummy          <- data.frame(ode(y = init, times = t, func = SEIR_model , parms = Opt_parDummy))
  fit_S.v[i,]       <- fitDummy[ , 2]
  fit_E.v[i,]       <- fitDummy[ , 3]
  fit_I_symp.v[i,]  <- fitDummy[ , 4]
  fit_I_asymp.v[i,] <- fitDummy[ , 5]
  Fit_I.v[i,]       <- fitDummy[ , 4] + fitDummy[ , 5]
  fit_cum_inf.v[i,] <- N - fitDummy[ , 2]
  
  fitted_incidence.v[i,]  <- p_symp_use *  fitDummy[ , 3] * eta
  #fitted_incidence.v[i,] <- beta(fitDummy[,1], delta = Opt_parDummy[1], epsilon = Opt_parDummy[2], theta = Opt_parDummy[3]) * fitDummy[ , 2] *  fitDummy[ , 4]/N 
  effective_reprod.v[i,] <- Basic_repr(fitDummy[,1], delta = Opt_parDummy[1], epsilon = Opt_parDummy[2], theta = Opt_parDummy[3], gamma = gammaD) * fitDummy[ , 2] /N
}


cum_inf_mean       <- apply(fit_cum_inf.v, MARGIN = 2, FUN = mean) 
cum_inf_median     <- apply(fit_cum_inf.v, MARGIN = 2, FUN = median) 
cum_inf_95_up_CRI  <- apply(fit_cum_inf.v, MARGIN = 2, FUN = CRI_95_up) 
cum_inf_95_low_CRI <- apply(fit_cum_inf.v, MARGIN = 2, FUN = CRI_95_low) 


fit_I_mean       <- apply(Fit_I.v, MARGIN = 2, FUN = mean) 
fit_I_median     <- apply(Fit_I.v, MARGIN = 2, FUN = median) 
fit_I_95_up_CRI  <- apply(Fit_I.v, MARGIN = 2, FUN = CRI_95_up) 
fit_I_95_low_CRI <- apply(Fit_I.v, MARGIN = 2, FUN = CRI_95_low) 

fitted_Incidence_mean       <- apply(fitted_incidence.v, MARGIN = 2, FUN = mean) 
fitted_Incidence_median     <- apply(fitted_incidence.v, MARGIN = 2, FUN = median) 
fitted_Incidence_95_up_CRI  <- apply(fitted_incidence.v, MARGIN = 2, FUN = CRI_95_up) 
fitted_Incidence_95_low_CRI <- apply(fitted_incidence.v, MARGIN = 2, FUN = CRI_95_low) 

effective_reprod_mean       <- apply(effective_reprod.v, MARGIN = 2, FUN = mean) 
effective_reprod_median     <- apply(effective_reprod.v, MARGIN = 2, FUN = median) 
effective_reprod_95_up_CRI  <- apply(effective_reprod.v, MARGIN = 2, FUN = CRI_95_up) 
effective_reprod_95_low_CRI <- apply(effective_reprod.v, MARGIN = 2, FUN = CRI_95_low) 




# Cummulative number of infected until 2020-04-11 and until 2020-05-01, with their 95% CI
# 2020-04-11 = dag 102
# 2020-05-01 = dag 122

as.Date(102, origin = "2019-12-31")
as.Date(122, origin = "2019-12-31")


fit_cum_inf[which(t==102)]

cum_inf_95_low_CRI[which(t==102)]
cum_inf_95_up_CRI[which(t==102)]

fit_cum_inf[which(t==102)]/N
cum_inf_95_low_CRI[which(t==102)]/N
cum_inf_95_up_CRI[which(t==102)]/N

fit_cum_inf[which(t==122)]
cum_inf_95_low_CRI[which(t==122)]
cum_inf_95_up_CRI[which(t==122)]

fit_cum_inf[which(t==122)]/N
cum_inf_95_low_CRI[which(t==122)]/N
cum_inf_95_up_CRI[which(t==122)]/N


as.Date(t[which(fit_I==max(fit_I))], origin = "2019-12-31" )
maxdagen <- as.Date(t[which(fit_I==max(fit_I))],origin = "2019-12-31" )

minDag <- as.Date(t[which(fit_I_95_low_CRI==max(fit_I_95_low_CRI))], origin = "2019-12-31" )
maxDag <- as.Date(t[which(fit_I_95_up_CRI==max(fit_I_95_up_CRI))], origin = "2019-12-31" )

as.Date(minDag, origin = "2019-12-31")
as.Date(maxDag, origin = "2019-12-31")




#############################################
## save estimated R0 and their uncertainty ##
#############################################

res_R0        <- round(c(p_asymp_use, p_lower_inf_use, R0_Mean[1], R0_Mean[length(Day)],effective_reprod_mean[length(Day)]), digits = 3)
names(res_R0) <- c("p_0", "q_0", "R0(start)", "R0(end)", "Re(end)")
CI_R01        <- paste("[",round(R0_low[1],digits=3), ", ", round(R0_high[1],digits=3),"]",sep="")
CI_ROend      <- paste("[",round(R0_low[2],digits=3), ", ", round(R0_high[2],digits=3),"]",sep="")
CI_Reend      <- paste("[",round(effective_reprod_95_low_CRI[length(Day)],digits=3), ", ", round(effective_reprod_95_up_CRI[length(Day)],digits=3),"]",sep="")
CIR0          <- c("", "", CI_R01, CI_ROend,CI_Reend)

MAT_R0    <- matrix(c(res_R0,CIR0), ncol =5, nrow=2, byrow=TRUE)
df.resR0  <- as.data.frame(MAT_R0)
colnames(df.resR0) <- names(res_R0)

XL_R0_file_name <- paste(table.path,"/Res_R0_p_non-reported_", p_asymp_use, "_infect_", p_lower_inf_use, ".xlsx", sep ="")
write.xlsx(df.resR0, XL_R0_file_name )



######################################
## save estimated days and their CI ##
######################################




fit_cum_low   <- cum_inf_95_low_CRI
fit_cum_high  <- cum_inf_95_up_CRI

res_days <-c(p_asymp_use, p_lower_inf_use, round(fit_cum_inf[which(t == 102)], digits=0), round(fit_cum_inf[which(t == 102)] / N, digits = 3) ,  
             round(fit_cum_inf[which(t == 122)],digits = 0), round(fit_cum_inf[which(t == 122)] / N, digits = 3))

res_maxprevalens <- round(max(fit_I), digits = 0)
res_days <- c(res_days, as.character(maxdagen), res_maxprevalens)

names(res_days)   <- c("p_0","q_0","2020-04-11", "", "2020-05-01","", "Peakdag", "Prevalens peakdag")
CIdag11mars       <- paste("[",round(fit_cum_low[which(t == 102)], digits=0), ", ", round(fit_cum_high[which(t == 102)], digits = 0), "]", sep="")
CIdag11marsProc   <- paste("[",round(fit_cum_low[which(t == 102)] / N, digits = 3), ", ",round(fit_cum_high[which(t == 102)]/N,digits=3),"]", sep="")
CIdag1maj         <- paste("[",round(fit_cum_low[which(t == 122)], digits=0), ", ", round(fit_cum_high[which(t == 122)],digits=0),"]", sep="")
CIdag1majProc     <- paste("[",round(fit_cum_low[which(t == 122)] / N, digits = 3),", ", round(fit_cum_high[which(t == 122)] / N , digits=3),"]", sep="")
CImaxdag          <- paste("[",as.character(as.Date(minDag, origin = "2019-12-31")), ", ", as.character(as.Date(maxDag, origin = "2019-12-31")), "]", sep="")
CIprevalensMaxdag <- paste("[",round(max(fit_I_95_low_CRI),digits = 0), ", ", round(max(fit_I_95_up_CRI), digits = 0), "]", sep="")



  
CI_dag  <- c("","", CIdag11mars, CIdag11marsProc, CIdag1maj, CIdag1majProc, CImaxdag, CIprevalensMaxdag)
MAT_dag <- matrix(c(res_days, CI_dag), ncol = 8, nrow = 2, byrow = TRUE)
df.dag  <- as.data.frame(MAT_dag)

colnames(df.dag) <- names(res_days)


XL_file_name <- paste(table.path,"/Res_dagar_p_non-reported_", p_asymp_use, "_infect_", p_lower_inf_use, ".xlsx", sep ="")
write.xlsx(df.dag, XL_file_name)


#----------------------------------------------------
# Save figures
#----------------------------------------------------


t   <- (Day[1]):(Day[length(Day)]+14) # time in days
fit <- data.frame(ode(y = init, times = t, func = SEIR_model, parms = Opt_par))

fit_S       <- fit[ , 2]
fit_E       <- fit[ , 3]
fit_I_symp  <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_I       <- fit_I_symp + fit_I_asymp
fit_I_E     <- fit_E + fit_I
fit_cum_inf <- N - fit_S



#fitted_incidence  <- beta(fit[,1], delta = Opt_par[1], epsilon = Opt_par[2], theta = Opt_par[3]) * fit_S *  fit_I_symp/N 
fitted_incidence  <- p_symp_use * fit_E * eta



fitted_incidence_low <- fitted_Incidence_95_low_CRI[1:length(t)]
fitted_incidence_high <- fitted_Incidence_95_up_CRI[1:length(t)]

fitted_I_high <- fit_I_95_up_CRI[1:length(t)]
fitted_I_low <- fit_I_95_low_CRI[1:length(t)]






NameNumber <- paste("/Incidence_number_infected_14Days_CI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".pdf",sep="")



pdf(paste(figure.path, NameNumber, sep=""), width=9*1.5, height=4.5*1.3 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 

## Prognosis until May



# plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
#      xlab="",ylim=c(0,roundUp(max(fitted_incidence))),xaxt='n')

plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(fitted_incidence_high))),xaxt='n')


grid(nx=NA, ny=NULL)


# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_april <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14)
NameDateMarchApril <- as.Date(dayatyear_march_april,origin ="2019-12-31")

abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

lines(fit$time, fitted_incidence_low, lty=2)
lines(fit$time, fitted_incidence_high, lty=2)


axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


## next plot

plot(fit$time, fit[ , 4] + fit[ , 5], ylim = c(0, roundUp(max(fitted_I_high)) ) , type="l", col="red", 
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)

lines(fit$time, fitted_I_high, lty=2)
lines(fit$time, fitted_I_low, lty=2)

axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)

dev.off()







NameNumberPNG <- paste("/Incidence_number_infected_14Days_CI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".png",sep="")

#width=9*1.5, height=4.5*1.3 

png(paste(figure.path, NameNumberPNG, sep=""),width = 9*1.5*75, height =  4.5*1.3*75 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 


plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(fitted_incidence_high))),xaxt='n')


grid(nx=NA, ny=NULL)


# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_april <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14)
NameDateMarchApril <- as.Date(dayatyear_march_april,origin ="2019-12-31")

abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

lines(fit$time, fitted_incidence_low, lty=2)
lines(fit$time, fitted_incidence_high, lty=2)



axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)



plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fitted_I_high)) ) , type="l", col="red", 
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)

lines(fit$time, fitted_I_high, lty=2)
lines(fit$time, fitted_I_low, lty=2)



axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)


dev.off()





NameNumber <- paste("/Incidence_number_infected_14Days_noCI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".pdf",sep="")

pdf(paste(figure.path, NameNumber, sep=""), width=9*1.5, height=4.5*1.3 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 

## Prognosis until May



plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(c(fitted_incidence,Observed_incidence)))),xaxt='n')



grid(nx=NA, ny=NULL)


# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_april <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14)
NameDateMarchApril <- as.Date(dayatyear_march_april,origin ="2019-12-31")

abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


## next plot

plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fit[ , 4]+fit[ , 5])) ) , type="l", col="red",
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)


axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)

dev.off()




NameNumberPNG <- paste("/Incidence_number_infected_14Days_noCI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".png",sep="")

#width=9*1.5, height=4.5*1.3 

png(paste(figure.path, NameNumberPNG, sep=""),width = 9*1.5*75, height =  4.5*1.3*75 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 


plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(c(fitted_incidence,Observed_incidence)))),xaxt='n')

grid(nx=NA, ny=NULL)

# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_april <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14)
NameDateMarchApril <- as.Date(dayatyear_march_april,origin ="2019-12-31")

abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)


axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)

## next plot

plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fit[ , 4]+fit[ , 5])) ) , type="l", col="red",
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')



grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)


axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)


dev.off()









#----------------------------------------------------
# Save figures 31 day forecast
#----------------------------------------------------


t   <- (Day[1]):(Day[length(Day)] + 31) # time in days
fit <- data.frame(ode(y = init, times = t, func = SEIR_model, parms = Opt_par))

fit_S       <- fit[ , 2]
fit_E       <- fit[ , 3]
fit_I_symp  <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_I       <- fit_I_symp + fit_I_asymp
fit_I_E     <- fit_E + fit_I
fit_cum_inf <- N - fit_S



#fitted_incidence  <- beta(fit[,1], delta = Opt_par[1], epsilon = Opt_par[2], theta = Opt_par[3]) * fit_S *  fit_I_symp/N 
fitted_incidence  <- p_symp_use * fit_E * eta

fitted_incidence_low <- fitted_Incidence_95_low_CRI[1:length(t)]
fitted_incidence_high <- fitted_Incidence_95_up_CRI[1:length(t)]

fitted_I_high <- fit_I_95_up_CRI[1:length(t)]
fitted_I_low <- fit_I_95_low_CRI[1:length(t)]





# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_may <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14 , 61 + 14 + 17 + 14 + 16, 61 + 14 + 17 + 14 + 16 + 14)
NameDateMarchMay <- as.Date(dayatyear_march_may,origin ="2019-12-31")


NameNumber <- paste("/Incidence_number_infected_31Days_CI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".pdf",sep="")



pdf(paste(figure.path, NameNumber, sep=""), width=9*1.5, height=4.5*1.3 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 

## Prognosis until May



# plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
#      xlab="",ylim=c(0,roundUp(max(fitted_incidence))),xaxt='n')

plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(fitted_incidence_high))),xaxt='n')


grid(nx=NA, ny=NULL)


abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

lines(fit$time, fitted_incidence_low, lty=2)
lines(fit$time, fitted_incidence_high, lty=2)


axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


## next plot

plot(fit$time, fit[ , 4] + fit[ , 5], ylim = c(0, roundUp(max(fitted_I_high)) ) , type="l", col="red", 
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)

lines(fit$time, fitted_I_high, lty=2)
lines(fit$time, fitted_I_low, lty=2)

axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)

dev.off()







NameNumberPNG <- paste("/Incidence_number_infected_31Days_CI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".png",sep="")

#width=9*1.5, height=4.5*1.3 

png(paste(figure.path, NameNumberPNG, sep=""),width = 9*1.5*75, height =  4.5*1.3*75 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 


plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(fitted_incidence_high))),xaxt='n')


grid(nx=NA, ny=NULL)



abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

lines(fit$time, fitted_incidence_low, lty=2)
lines(fit$time, fitted_incidence_high, lty=2)



axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)



plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fitted_I_high)) ) , type="l", col="red", 
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)

lines(fit$time, fitted_I_high, lty=2)
lines(fit$time, fitted_I_low, lty=2)



axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)


dev.off()





NameNumber <- paste("/Incidence_number_infected_31Days_noCI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".pdf",sep="")

pdf(paste(figure.path, NameNumber, sep=""), width=9*1.5, height=4.5*1.3 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 

## Prognosis until May



plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(c(fitted_incidence,Observed_incidence)))),xaxt='n')



grid(nx=NA, ny=NULL)


abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


## next plot

plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fit[ , 4]+fit[ , 5])) ) , type="l", col="red",
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)


axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)

dev.off()




NameNumberPNG <- paste("/Incidence_number_infected_31Days_noCI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".png",sep="")

#width=9*1.5, height=4.5*1.3 

png(paste(figure.path, NameNumberPNG, sep=""),width = 9*1.5*75, height =  4.5*1.3*75 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 


plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(c(fitted_incidence,Observed_incidence)))),xaxt='n')

grid(nx=NA, ny=NULL)


abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)


axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)

## next plot

plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fit[ , 4]+fit[ , 5])) ) , type="l", col="red",
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')



grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)


axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)


dev.off()








