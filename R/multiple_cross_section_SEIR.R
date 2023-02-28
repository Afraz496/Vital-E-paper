# change your working directory here 
HOME_WD <- 'C:/Users/CH3079_SA/Desktop/Vital-E/Phase 1'
data_dir <- '/virosolver_paper/Vital-E'
filename <- '/control.csv'
library(virosolver)
library(tidyverse)
library(patchwork)
library(lazymcmc)
library(foreach)
library(doParallel)

# source plotting functions
source('plottingfuncs.R')

# Parallelization
cl <- parallel::makeCluster(4, setup_strategy = "sequential")
registerDoParallel(cl)

############### DATA #################
##############*******#################

# Load in the example file
data(example_seir_partab)

# Load in the Dental Data
data <- read_csv(paste0(HOME_WD, data_dir, filename))

dataRAW <- read_csv(paste0(HOME_WD, data_dir, filename))

# Clean the data with missing values
# mutate missing values, and modify the dataframe
# if Ct Value is missing -> 40, if it is 0 -> 40
# data <- data %>%
#   mutate(Ct.e = replace(Ct.e,is.na(Ct.e),40)) %>%
#   mutate(Ct.e = replace(Ct.e, Ct.e==0, 40))
# 
# dataRAW <- dataRAW %>%
#   mutate(Ct.e = replace(Ct.e,is.na(Ct.e),40)) %>%
#   mutate(Ct.e = replace(Ct.e, Ct.e==0, 40))

# Fix column names for MCMC framework

data <- data %>%
  mutate(collection_date = as.Date(collection_date, "%m/%d/%Y"))

dataRAW <- dataRAW %>%
  mutate(collection_date = as.Date(collection_date, "%m/%d/%Y"))

data <- data %>% 
  filter(collection_date > as.Date('2020-10-01'))

dataRAW <- dataRAW %>% 
  filter(collection_date > as.Date('2020-10-01'))

# Make Calendar Time
data <- data %>%
  arrange(collection_date) %>%
  mutate(first_date = min(collection_date, na.rm=TRUE)) %>%
  mutate(calendar_time = collection_date - first_date) %>%
  mutate(calendar_time = as.numeric(calendar_time))

dataRAW <- dataRAW %>%
  arrange(collection_date) %>%
  mutate(first_date = min(collection_date, na.rm=TRUE)) %>%
  mutate(calendar_time = collection_date - first_date) %>%
  mutate(calendar_time = as.numeric(calendar_time))

dataRAW <- dataRAW %>%
  mutate(t = calendar_time,
         ct = rnasep_correction) 

data <- data %>%
  mutate(t = calendar_time,
         ct = rnasep_correction) 

# Prepare columns for pipeline
data <- data %>%
  select(t, ct)

posdataRAW <- dataRAW %>%
  filter(!is.na(ct))

# Select Horizon Here
data <- data %>%
  filter(t %in% c(61,69,75)) %>% 
  filter(!is.na(ct))

#----- VITAL-E DATA MODIFICATIONS-----#
## Look at the code to see how the prior is getting sample over the prior
example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='t0', values = 0))

example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='t0', fixed = 1))

example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='t0', lower_bound = -10))

example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='I0', fixed = 0))

example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='I0', values = 3.20512821e-6))

example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='infectious', fixed = 1))

example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='incubation', fixed = 1))

example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='I0', lower_bound = 6e-9))

example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='I0', upper_bound = 0.1))

example_seir_partab <- example_seir_partab %>%
  rows_update(tibble(names='viral_peak', values = 15))


## Fix everything but I0 and R0 
# Vary the rest 1 by 1 to see what kind of SEIR curve we get
######################################

# FOR THIS MODEL SHIFT ALL THE CALENDAR DAYS by an amount and filter to > 0
t_shift <- 54 # we need our first horizon to be on day 7
data <- data %>% 
  mutate(t = t - t_shift) %>% 
  filter(t > 0)

dataRAW <- dataRAW %>% 
  mutate(t = t-t_shift) %>% 
  filter(t > 0)
## ---Multiple Cross Section Model (Vignette 4.5) ---#
## MCMC chain options
mcmc_pars <- c("iterations"=500000,"popt"=0.44,"opt_freq"=2000,
               "thin"=2500,"adaptive_period"=200000,"save_block"=100)

## Set pointer to the SEIR model as the incidence function
incidence_function <- solveSEIRModel_lsoda_wrapper

## This is for the GP version
times <- 0:max(data$t)
mat <- matrix(rep(times, each=length(times)),ncol=length(times))
t_dist <- abs(apply(mat, 2, function(x) x-times))
par_tab <- example_seir_partab

pars <- par_tab$values
names(pars) <- par_tab$names

## Pull out the current values for each parameter, and set these as the prior means
means <- par_tab$values
names(means) <- par_tab$names

## Set standard deviations of prior distribution
sds_seir <- c("obs_sd"=0.5,"viral_peak"=2,
            "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
            "prob_detect"=0.03,
            "incubation"=0.25, "infectious"=0.5)

## Define a function that returns the log prior probability for a given vector of parameter
## values in `pars`, given the prior means and standard deviations set above.
prior_func_seir <- function(pars,...){
  ## Ct model priors
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_seir["obs_sd"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_seir["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_seir["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_seir["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_seir["level_switch"],log=TRUE)
  ## Beta prior on the prob_detect parameter to ensure between 0 and 1
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_seir["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  ## SEIR model priors
  incu_prior <- dlnorm(pars["incubation"],log(means[which(names(means) == "incubation")]), sds_seir["incubation"], TRUE)
  infectious_prior <- dlnorm(pars["infectious"],log(means[which(names(means) == "infectious")]),sds_seir["infectious"],TRUE)
  
  ## Sum up
  obs_sd_prior + viral_peak_prior + 
    wane_2_prior + tswitch_prior + level_prior + beta_prior +
    incu_prior + infectious_prior
}

posterior_function <- create_posterior_func(parTab=par_tab, 
                                            data=data, 
                                            PRIOR_FUNC=prior_func_seir, 
                                            INCIDENCE_FUNC=incidence_function,
                                            t_dist=t_dist, use_pos = FALSE)
posterior_function(par_tab$values)


dir.create("mcmc_chains/correction/readme_multiple_cross_section",recursive=TRUE)
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
                                          prior_func_seir)
  
  output <- run_MCMC(parTab=start_tab,
                     data=data,
                     INCIDENCE_FUNC=incidence_function,
                     PRIOR_FUNC=prior_func_seir,
                     mcmcPars=mcmc_pars,
                     filename=paste0("mcmc_chains/correction/readme_multiple_cross_section/readme_seir_",chain_no),
                     CREATE_POSTERIOR_FUNC=create_posterior_func,
                     mvrPars=NULL,
                     OPT_TUNING=0.2,
                     use_pos=FALSE,
                     t_dist=t_dist)
}

cts <- rep(virosolver::viral_load_func(pars, seq(0,60,0.1), FALSE, 0),1)

ct_line_plot <- tibble(cts = cts, t = seq(0,60,0.1)) %>%
  ggplot() + 
  geom_line(aes(x=t,y=cts)) + 
  theme_bw() +
  xlab("Time since infected") +
  ylab("Ct")
ct_line_plot

violin <- ggplot(data) +
  geom_violin(aes(x=as.factor(t),y=ct), draw_quantiles = c(0.5)) +
  xlab('Time') +
  ylab('Ct')

violin
## Read in the MCMC chains
chains <- load_mcmc_chains(location="mcmc_chains/correction/readme_multiple_cross_section",
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
chains_melted <- chains$chain %>% as_tibble %>% group_by(chain) %>% mutate(sampno=1:n()) %>% pivot_longer(-c(sampno,chain))
## Look at trace plots
p_trace_gp <- chains_melted %>%
  filter(!(name %in% paste0("prob.",1:max(times)))) %>%
  ggplot() + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~name,scales="free_y") + 
  scale_color_viridis_d(name="Chain") + 
  theme_bw() +
  xlab("Iteration") +
  ylab("Value")
p_trace_gp

## Load in MCMC chains again, but this time read in the fixed parameters too 
## to ensure that the posterior draws are compatible with the model functions
chains <- load_mcmc_chains(location="mcmc_chains/correction/readme_multiple_cross_section/",
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
## Violin Plots for all the priors
violin <- plot_violins_mcmc_params(chains_melted)

# add the priors to the plot
violin

chain_comb <- as.data.frame(chains$chain)
chain_comb$sampno <- 1:nrow(chain_comb)



## Load in true incidence curve to compare to our prediction
data(example_seir_incidence)
predictions <- plot_prob_infection(chain_comb,nsamps=10, INCIDENCE_FUNC=incidence_function,
                                   solve_times=0:max(data$t),obs_dat=data,
                                   true_prob_infection=example_seir_incidence,
                                   smooth=TRUE) ## Smooth the trajectories a bit
p_incidence_prediction <- predictions$plot + scale_x_continuous(limits=c(0,40))
p_incidence_prediction

# plot incidence vs prediction (NOT PROBABILITY)

# Load in Vital-E incidence curve
# reconstruct incidence from linelist
get_incidence_from_linelist <- function(data){
  data %>%
    group_by(t) %>%
    summarise(incidence = n())
}


incidence <- get_incidence_from_linelist(dataRAW)

incidence <- incidence %>%
  mutate(prob_infection = incidence,
         t = t)

posdataRAW <- dataRAW %>% 
  filter(!is.na(rnasep_correction))

dataRAW <- dataRAW %>%
  filter(t < 50)

posdataRAW <- posdataRAW %>% 
  filter(t < 50)
source('plottingfuncs_anon.R')
start_date <- min(dataRAW$collection_date)
predictions <- plot_infection_incidence_anonymous(chain_comb, 
                                        nsamps=100,
                                        posdata = posdataRAW, rawdata = dataRAW,
                                        INCIDENCE_FUNC=incidence_function,
                                        solve_times=0:max(data$t),
                                        obs_dat=data,
                                        true_prob_infection=incidence, population = nrow(posdataRAW),
                                        horizons=FALSE,
                                        mode='Incidence', start_date = start_date)
predictions


## Use create_posterior_func to return the predicted Ct distribution rather than the posterior probability
model_func_gp <- create_posterior_func(par_tab,data,NULL,incidence_function,"model")
## Pass model_func to a plotting function to observe predicted Ct distribution against data
p_distribution_fit_gp <- plot_distribution_fits(chain_comb, data, model_func_gp,100,pos_only=FALSE)
## Joining, by = "t"
## Joining, by = c("t", "sampno")
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
p_distribution_fit_gp

