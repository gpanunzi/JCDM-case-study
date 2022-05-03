
# Packages ----------------------------------------------------------------

require(tidyverse)
require(magrittr)
require(lubridate)
require(RMark)
require(R2ucare)
require(zeallot)

# Loading Data --------------------------------------------------------------------

# Sight Matrix
load("Data/CapDolph_SightMatrix_2017_2020.RData")

# Effort by days
load("Data/CapDolph_EffortByDays.RData")

# Clustered data
load("Data/SightMatrix_Clustered.RData")


# Manipulating data -------------------------------------------------------

# Effort days per year
effDays <- effort_days %>% unite(YM, Y, M, sep="-") %>% 
  mutate(YM=as.Date(paste(YM,"-01",sep=""))) 
effort_ndays <- effort_days %$% n

# Year-Months date
YM <- raw_data_turs %>% group_by(year(Data), month(Data)) %>% 
  summarise(fDat = first(Data)) %>% mutate(fDat=floor_date(fDat, "month")) %$% fDat

# Number of months per year
MpY <- raw_data_turs %>%
  mutate(Y=year(Data), M=month(Data)) %>% distinct(Y, M) %>% group_by(Y) %>% 
  summarise(n=n()) %$% n

# Residency groups: 1=resident, 2=non resident, 3=transient
ClusByID <- data_turs_clustered %>% distinct(ID, Group) %>% mutate(Group=as.factor(Group))

# Deriving the sighting matrix with group labels
mat_avv <- raw_data_turs %>% 
  filter(`Mark level`=="WM") %>% 
  mutate(YM=floor_date(Data, "month")) %>% group_by(ID, YM) %>% 
  summarise(Observed=as.numeric(sum(Observed)>0)) %>% 
  spread(., YM, Observed) %>% 
  left_join(ClusByID) %>% 
  unite("ch", !c("ID", "Group"), sep="", remove=F) %>% 
  column_to_rownames("ID")


# Time_intervals between months (i.e. occasions)
t_list <- c(2,1,1,         # 2017
            6,1,1,1,       # 2018
            7,1,1,1,1,1,1, # 2019
            9,1,1,1,1      # 2020
)

# Quantities of interest ---------------------------------------------------

# Total number of sighted animals
n <- raw_data_turs %>% filter(Observed==1) %>% 
  distinct(ID, `Mark level`) %>% 
  summarise(n=n()) %$% n

# Total number of sighted animals per year and month
nYM <- raw_data_turs %>% filter(Observed==1) %>% 
  group_by(year(Data), month(Data)) %>% distinct(ID, `Mark level`) %>% 
  summarise(n=n()) %$% n

# Total number of sighted animals per year
nY <- raw_data_turs %>% filter(Observed==1) %>% 
  group_by(Y=year(Data)) %>% distinct(ID, Y) %>% 
  summarise(n=n()) %$% n %>% rep(., times=MpY)

# Well Marked rate
theta <- raw_data_turs %>% filter(Observed==1) %>% 
  distinct(ID, `Mark level`) %>% 
  filter(`Mark level`=="WM") %>% 
  summarise(n=n()) %$% n/n

# Well marked rate per year and month
thetaYM <- raw_data_turs %>% filter(Observed==1) %>% 
  group_by(year(Data), month(Data)) %>% distinct(ID, `Mark level`) %>% 
  filter(`Mark level`=="WM") %>% 
  summarise(n=n()) %$% n/nYM

# Well marked rate per year
thetaY <- raw_data_turs %>% filter(Observed==1) %>% 
  group_by(Y=year(Data)) %>% distinct(ID, `Mark level`) %>% 
  filter(`Mark level`=="WM") %>% 
  summarise(n=n()) %$% n %>% rep(., times=MpY)/nY


# C Hat - Closure test -------------------------------------------------------------------

# CJS test
X_model <- mat_avv %>% select(-Group, -ch) %>%  as.matrix()
overall_cjs <-  overall_CJS(X = X_model, freq = rep(1, nrow(X_model)))
overall_cjs

# Correction by degrees of freedom
chi <- overall_cjs["chi2"] %>% as.numeric()
dof <- overall_cjs["degree_of_freedom"] %>% as.numeric()
(c_hat <- chi/dof)


# POPAN fit -------------------------------------------------------------------

# Processing data for POPAN
dolph <- process.data(data = as.data.frame(mat_avv %>% select(ch, Group)), 
                      model = "POPAN", groups="Group",
                      time.intervals = t_list, begin.time = 1)
dolph.ddl <- make.design.data(dolph)

# Setting up the effort as a covariate for the detection
dolph.ddl$p$effort <- effort_ndays

# Auxiliary function for running different models
run_models <- function(dolph2, dolph2.ddl)
{
  # Detection prob formulas
  p.ct <- list(formula=~1)
  p.logeff <- list(formula=~log(effort))
  p.time <- list(formula=~time)
  plist <- list(p=p.ct, p=p.logeff, p=p.time)
  
  # Survival prob formulas
  phi.ct <- list(formula=~1)
  phi.time <- list(formula=~time)
  phi.group <- list(formula=~group)
  philist <- list(phi=phi.ct, phi=phi.time, phi=phi.group)

  # Entrance prob formulas  
  pent.ct <- list(formula=~1)
  pent.time <- list(formula=~time)
  pent.group <- list(formula=~group)
  pentlist <- list(pent=pent.ct, pent=pent.time, pent=pent.group)
  
  # Building all combinations
  combs <- cross3(1:length(plist), 1:length(philist), 1:length(pentlist)) %>% 
    reduce(., rbind)
  
  # Fitting the various models through mark
  mod <- list()
  for (i in 1:nrow(combs))
  {
    idx <- combs[i, ] %>% reduce(., c)
    modIdx <- paste("Mod", idx[1], idx[2], idx[3], sep="", collapse = "")
    mod[[modIdx]] <- mark(dolph2, dolph2.ddl, model.parameters=list(Phi=philist[[idx[2]]], 
                                                                    p=plist[[idx[1]]], 
                                                                    pent=pentlist[[idx[3]]]), 
                          output = F, delete=T, silent=T)
  }
  
  # Saving all the output in the global environment
  list2env(mod, .GlobalEnv)
}


# Executing the function to run the various models (it may take some minutes, ignore the warnings)
run_models(dolph, dolph.ddl)

# Collecting all the models and comparing the results
modColl <- collect.models()
adjust.chat(chat = c_hat, modColl)

# Selecting the output from the best model
bestMod <- modColl$Mod233
bestModResults <- popan.derived(dolph, bestMod, normal = F)

# Results from best model -----------------------------------------------------------------

# Estimated coefficients
(bestCoeffs <- bestMod$results$beta)

# Detection, enrance and suvival for different groups at different occasions
(pEst <- bestMod$results$real$estimate[grep("^p ", row.names(bestMod$results$real))])
(pentEst <- bestMod$results$real$estimate[grep("^pent ", row.names(bestMod$results$real))])
(phiEst <- bestMod$results$real$estimate[grep("^Phi ", row.names(bestMod$results$real))])

# Gross estimate of the population size with intervals

# Residents
(N_estRes <- bestModResults$NGross[1])
(N_seRes <- sqrt(bestModResults$NGross.vcv[1,1]))
alpha <- 0.05
cRes <- exp(qnorm(1-alpha/2)*sqrt(log(1+(N_seRes^2)/(N_estRes^2))))

# Not residents
(N_estNRes <- bestModResults$NGross[2])
(N_seNRes <- sqrt(bestModResults$NGross.vcv[2,2]))
alpha <- 0.05
cNRes <- exp(qnorm(1-alpha/2)*sqrt(log(1+(N_seNRes^2)/(N_estNRes^2))))

# Part-time residents
(N_estPRes <- bestModResults$NGross[3])
(N_sePRes <- sqrt(bestModResults$NGross.vcv[3,3]))
cPRes <- exp(qnorm(1-alpha/2)*sqrt(log(1+(N_sePRes^2)/(N_estPRes^2))))

# Total population (standard error includes covariance)
(N_est <- N_estRes+N_estNRes+N_estPRes)
(N_se <- sqrt(bestModResults$NGross.vcv[1,1]+bestModResults$NGross.vcv[2,2]+bestModResults$NGross.vcv[3,3]+
                2*bestModResults$NGross.vcv[1,2]+
                2*bestModResults$NGross.vcv[1,3]+
                2*bestModResults$NGross.vcv[2,3]))
c_ <- exp(qnorm(1-alpha/2)*sqrt(log(1+(N_se^2)/(N_est^2))))

# Final table with all the point estimates by group (and total)
data.frame("Lower" = c(N_est/c_, N_estRes/(cRes), N_estNRes/(cNRes), N_estPRes/(cPRes)), 
           "Ntot" = c(N_est, N_estRes, N_estNRes, N_estPRes), 
           "Upper" = c(N_est*c_, N_estRes*(cRes), N_estNRes*(cNRes), N_estPRes*(cPRes))) %>% 
  round %>% set_rownames(c("Total", "Residents", "Transient", "Part-Time"))

# Total population by occasion
(N_est_ByOcc <- bestModResults$Nbyocc$N)
(N_se_ByOcc <- bestModResults$Nbyocc$se)
# Plot
bestModResults$Nbyocc %>% mutate(date = YM) %>% ggplot() + geom_point(aes(x=date, y=N)) + 
  geom_segment(aes(x=date, xend=date, y=LCL, yend=UCL)) +
  theme_light()


# Reporting to total by accounting for the mark rate ---------------------------------------------------------

# Auxiliary function to correct the estimates for the mark rate (Delta Method)
compute_Ntot <- function(estimate, stderr, theta, N_observed, alpha = 0.05){
  
  zalpha <- qnorm(1 - alpha/2)
  
  Ntot <- estimate/theta
  se_Ntot <- sqrt(Ntot^2*(((stderr)/estimate)^2+((1-theta)/(N_observed*theta))))
  c_hat_corr <- exp(zalpha*sqrt(log(1+(se_Ntot/Ntot)^2)))
  
  return(list(Ntot = Ntot, se_Ntot = se_Ntot, c_hat = c_hat_corr))
}

# To the universe Total
(c(Ntot_est, Ntot_se, c_hat) %<-% compute_Ntot(estimate = N_est, stderr = N_se, 
                                               theta = theta, N_observed = n))

# To the universe Resident
c(NtotRes_est, NtotRes_se, c_hatRes) %<-% compute_Ntot(estimate = N_estRes, stderr = N_seRes, 
                                                       theta = theta, N_observed = n)

# To the universe Non Resident
c(NtotNRes_est, NtotNRes_se, c_hatNRes) %<-% compute_Ntot(estimate = N_estNRes, stderr = N_seNRes, 
                                                          theta = theta, N_observed = n)

# Final table with all the point estimates by group (and total)
data.frame("Lower" = c(Ntot_est/(c_hat), NtotRes_est/(c_hatRes), NtotNRes_est/(c_hatNRes)), 
           "Ntot" = c(Ntot_est, NtotRes_est, NtotNRes_est), 
           "Upper" = c(Ntot_est*(c_hat), NtotRes_est*(c_hatRes), NtotNRes_est*(c_hatNRes))) %>% 
  round %>% set_rownames(c("Total", "Residents", "Non Residents"))

# Total population by occasion
(NTot_est_ByOcc <- compute_Ntot(estimate = N_est_ByOcc, stderr = N_se_ByOcc, 
                                theta = thetaY, N_observed = nY) %>% 
    reduce(., .f=bind_cols) %>% set_colnames(c("est", "se", "chat")) %>% 
    mutate(date = YM))
# Plot
NTot_est_ByOcc %>% ggplot() + geom_point(aes(x=date, y=est)) + 
  geom_segment(aes(x=date, xend=date, y=est/chat, yend=est*chat)) + 
  theme_light()

