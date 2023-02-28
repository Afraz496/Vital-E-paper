library(virosolver)
library(tidyverse)
library(patchwork)
library(lazymcmc)
library(foreach)
library(doParallel)

# change your working directory here
HOME_WD <- 'D:/Job/BCCDC/Vital-E/Phase 1'

# source plotting functions
source('plottingfuncs.R')

############### DATA #################
##############*******#################
# Load the Dental Conference Data
data_dir <- '/virosolver_paper/Vital-E'

## Read in the GP model parameter control table
data(example_gp_partab)

# Load in the Dental Data
data <- read_csv(paste0(HOME_WD, data_dir, '/tabor.csv'))

dataRAW <- read_csv(paste0(HOME_WD, data_dir, '/tabor.csv'))

# Clean the data with missing values
# mutate missing values, and modify the dataframe
# if Ct Value is missing -> 40, if it is 0 -> 40
data <- data %>%
  mutate(Ct.e = replace(Ct.e,is.na(Ct.e),40)) %>%
  mutate(Ct.e = replace(Ct.e, Ct.e==0, 40))

dataRAW <- dataRAW %>%
  mutate(Ct.e = replace(Ct.e,is.na(Ct.e),40)) %>%
  mutate(Ct.e = replace(Ct.e, Ct.e==0, 40))

# Fix column names for MCMC framework

dataRAW <- dataRAW %>%
  mutate(t = calendar_time,
         ct = Ct.e) 

data <- data %>%
  mutate(t_raw = calendar_time,
         ct = Ct.e) 

# Unique to Tabor home but use this if you have multiple outbreaks in your data
data <- data %>%
  filter(t_raw >= 0, t_raw <= 28) %>%
  mutate(t = t_raw) %>%
  select(t, ct)

# Filter for Tabor Home specifically
dataRAW <- dataRAW %>%
  filter(t >= 0, t <= 28)

posdataRAW <- dataRAW %>%
  filter(ct != 40)

# Load in Phase 2 data
phase2 <- read_csv(paste0(HOME_WD, data_dir, '/phase2.csv'))


phase2_RAW <- read_csv(paste0(HOME_WD, data_dir, '/phase2.csv'))

#----- VITAL-E DATA MODIFICATIONS-----#
example_gp_partab <- example_gp_partab %>%
  rows_update(tibble(names='t0', values = 0))

######################################

phase2 <- phase2 %>%
  mutate(t = calendar_time,
         ct = Ct.e) %>%
  filter(!is.na(ct),
         !is.na(t))

phase2_RAW <- phase2_RAW %>%
  mutate(t = calendar_time,
         ct = Ct.e) %>%
  filter(!is.na(ct),
         !is.na(t))

phase2 <- phase2 %>%
  select(t, ct) 

phase2 <- phase2 %>%
  filter(t %in% c(10, 22, 35, 46, 59, 66, 74, 80))

## Shift time in data (if required)
time_shift <- 2
dataRAW <- dataRAW %>%
  mutate(t = t+time_shift)

# Select the Time slices here
data <- dataRAW %>%
  select(t,ct) %>%
  filter(t %in% c(2, 4, 11, 9, 15, 16))
## ---Multiple Cross Section Model (Vignette 4.5) ---#
## MCMC chain options
mcmc_pars <- c("iterations"=500000,"popt"=0.44,"opt_freq"=2000,
               "thin"=2500,"adaptive_period"=200000,"save_block"=100)

## Set pointer to the Gaussian Process model as the incidence function
incidence_function <- gaussian_process_model


## This is for the GP version
times <- 0:max(data$t)
mat <- matrix(rep(times, each=length(times)),ncol=length(times))
t_dist <- abs(apply(mat, 2, function(x) x-times))
par_tab <- example_gp_partab
par_tab <- bind_rows(par_tab[par_tab$names != "prob",], par_tab[par_tab$names == "prob",][1:length(times),])
pars <- par_tab$values
names(pars) <- par_tab$names

## Pull out the current values for each parameter, and set these as the prior means
means <- par_tab$values
names(means) <- par_tab$names

## Set standard deviations of prior distribution
sds_gp <- c("obs_sd"=0.5,"viral_peak"=2,
            "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
            "prob_detect"=0.03,
            "incubation"=0.25, "infectious"=0.5)

## Define a function that returns the log prior probability for a given vector of parameter
## values in `pars`, given the prior means and standard deviations set above.
## Prior for GP version
prior_func_gp <- function(pars, ...){
  par_names <- names(pars)
  
  ## Viral kinetics parameters
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_gp["obs_sd"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_gp["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_gp["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_gp["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_gp["level_switch"],log=TRUE)
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_gp["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  ### VERY IMPORTANT
  ## Gaussian process prior, un-centered version
  k <- pars[which(par_names=="prob")]
  ## Leave this - correct for uncentered version as per Chapter 14 Statistical Rethinking
  prob_priors <- sum(dnorm(k, 0, 1, log=TRUE))
  #########
  
  nu_prior <- dexp(pars["nu"], 1/means[which(names(means) == "nu")],log=TRUE)
  rho_prior <- dexp(pars["rho"], 1/means[which(names(means) == "rho")],log=TRUE)
  
  obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior +
    level_prior + beta_prior + prob_priors +
    nu_prior + rho_prior
}

posterior_function <- create_posterior_func(parTab=par_tab, 
                                            data=data, 
                                            PRIOR_FUNC=prior_func_gp, 
                                            INCIDENCE_FUNC=incidence_function,
                                            t_dist=t_dist, use_pos = FALSE)
posterior_function(par_tab$values)


dir.create("mcmc_chains/readme_multiple_cross_section",recursive=TRUE)
##################################
## RUN THE MCMC FRAMEWORK
## Run 3 MCMC chains. Note that it is possible to parallelize this loop with foreach and doPar
## Note the `use_pos` argument needs to be set here too
nchains <- 3
res <- foreach(chain_no=1:nchains,.packages = c("virosolver","lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
  ## Get random starting values
  start_tab <- generate_viable_start_pars(par_tab,data,
                                          create_posterior_func,
                                          incidence_function,
                                          prior_func_gp)
  
  output <- run_MCMC(parTab=start_tab,
                     data=data,
                     INCIDENCE_FUNC=incidence_function,
                     PRIOR_FUNC=prior_func_gp,
                     mcmcPars=mcmc_pars,
                     filename=paste0("mcmc_chains/readme_multiple_cross_section/readme_gp_",chain_no),
                     CREATE_POSTERIOR_FUNC=create_posterior_func,
                     mvrPars=NULL,
                     OPT_TUNING=0.2,
                     use_pos=FALSE,
                     t_dist=t_dist)
}

## Prototype for Prior Function (before MCMC)
start_tab <- generate_viable_start_pars(par_tab,data,
                                        create_posterior_func,
                                        incidence_function,
                                        prior_func_gp)

# Treat this the same as the chains variable
before_mcmc_melted <- 

## Read in the MCMC chains
chains <- load_mcmc_chains(location="mcmc_chains/readme_multiple_cross_section",
                           parTab=par_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=TRUE,
                           unfixed=TRUE,
                           multi=FALSE)
## [1] "mcmc_chains/readme_multiple_cross_section/readme_gp_1_univariate_chain.csv"
## [2] "mcmc_chains/readme_multiple_cross_section/readme_gp_2_univariate_chain.csv"
## [3] "mcmc_chains/readme_multiple_cross_section/readme_gp_3_univariate_chain.csv"
## [[1]]
## [1] 201
## 
## [[2]]
## [1] 201
## 
## [[3]]
## [1] 201
## Reshape for plotting
chains_melted_RAW <- chains$chain %>% as_tibble %>% group_by(chain) %>% mutate(sampno=1:n()) %>% pivot_longer(-c(sampno,chain))
## Look at trace plots
p_trace_gp <- chains_melted_RAW %>%
  filter(!(name %in% paste0("prob.",1:max(times)))) %>%
  ggplot() + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~name,scales="free_y") + 
  scale_color_viridis_d(name="Chain") + 
  theme_bw() +
  xlab("Iteration") +
  ylab("Value")
p_trace_gp

# Look at Violin Plots
chains_melted <- chains_melted_RAW %>%
  filter(!(name %in% paste0("prob.",1:max(times))))
violin <- plot_violins_mcmc_params(chains_melted)
violin

## Load in MCMC chains again, but this time read in the fixed parameters too 
## to ensure that the posterior draws are compatible with the model functions
chains <- load_mcmc_chains(location="mcmc_chains/readme_multiple_cross_section/",
                           parTab=par_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=FALSE,
                           unfixed=FALSE,
                           multi=FALSE)
## [1] "mcmc_chains/readme_multiple_cross_section//readme_gp_1_univariate_chain.csv"
## [2] "mcmc_chains/readme_multiple_cross_section//readme_gp_2_univariate_chain.csv"
## [3] "mcmc_chains/readme_multiple_cross_section//readme_gp_3_univariate_chain.csv"
## [[1]]
## [1] 201
## 
## [[2]]
## [1] 201
## 
## [[3]]
## [1] 201
## Do some reshaping to allow correct subsampling (we need each sampno to correspond to one unique posterior draw)
chain_comb <- as.data.frame(chains$chain)
chain_comb$sampno <- 1:nrow(chain_comb)

## Load in true incidence curve to compare to our prediction
predictions <- plot_prob_infection(chain_comb,nsamps=100, INCIDENCE_FUNC=incidence_function,
                                   solve_times=0:max(data$t),obs_dat=data,
                                   true_prob_infection=NULL,
                                   smooth=TRUE) ## Smooth the trajectories a bit
p_incidence_prediction <- predictions$plot + scale_x_continuous(limits=c(0,28))
p_incidence_prediction


# plot incidence vs prediction (NOT PROBABILITY)

# Load in Vital-E incidence curve
# reconstruct incidence from linelist
get_incidence_from_linelist <- function(data){
  data %>%
    group_by(t) %>%
    summarise(incidence = n())
}


incidence <- get_incidence_from_linelist(data)

incidence <- incidence %>%
  mutate(prob_infection = incidence,
         t = t)


predictions <- plot_infection_incidence(chain_comb, 
                                        nsamps=100,
                                        posdata = posdataRAW,
                                        rawdata = dataRAW,
                                        INCIDENCE_FUNC=incidence_function,
                                        solve_times=0:max(data$t),
                                        obs_dat=NULL,
                                        true_prob_infection=incidence, population = NULL)

p_incidence_prediction <- predictions$plot + scale_x_continuous(limits=c(0,28)) + scale_y_continuous(limits=c(0,50))
p_incidence_prediction

## Use create_posterior_func to return the predicted Ct distribution rather than the posterior probability
model_func_gp <- create_posterior_func(par_tab,phase2,NULL,incidence_function,"model")
## Pass model_func to a plotting function to observe predicted Ct distribution against data
p_distribution_fit_gp <- plot_distribution_fits(chain_comb, phase2, model_func_gp,100,pos_only=FALSE)
## Joining, by = "t"
## Joining, by = c("t", "sampno")
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
p_distribution_fit_gp

