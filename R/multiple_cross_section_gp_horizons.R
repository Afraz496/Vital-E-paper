# change your working directory here
HOME_WD <- 'C:/Users/CH3079_SA/Desktop/Vital-E/Phase 1'

library(virosolver)
library(tidyverse)
library(patchwork)
library(lazymcmc)
library(foreach)
library(doParallel)
library(readr)
library(janitor)

filter_horizon_data <- function (data, t_shift, horizon){
  data <- data %>%
    mutate(t = t+t_shift)
  
  horizon <- horizon + t_shift
  
  data <- data %>%
    filter(t %in% horizon)
  
  return(data)
} 

# Parallelization
cl <- parallel::makeCluster(4, setup_strategy = "sequential")
registerDoParallel(cl)

# source plotting functions
source('plottingfuncs.R')

# Read the data once to decide horizons
data_dir <- '/virosolver_paper/Vital-E'
dataRAW <- read_csv(paste0(HOME_WD, data_dir, '/phase2_1_first_pos_msp.csv'), col_types=cols(Ct.e = "d", Ct.Rdrp = "d"))

# Set when you wanna start your data observations
beginning_times <- 258
end_times <- 350

# Change horizions here
# NOTE: These are horizons BEFORE shifting
#horizons <- list(seq(beginning_times,end_times,by=75), seq(beginning_times,end_times,by=50), seq(beginning_times,end_times,by=25))
horizons <- list(c(258, 274, 315), c(258, 274, 315, 336), c(258, 274, 315, 336, 350))
t_shift<-5
for (i in 1:length(horizons)){
  horizons[[i]]<- horizons[[i]]+t_shift
}

# Change the SEIR PARTAB
## Read in the GP model parameter control table
data(example_gp_partab)

example_gp_partab <- within(example_gp_partab, fixed[names == 'viral_peak'] <- 1)
example_gp_partab <- within(example_gp_partab, fixed[names == 'obs_sd'] <- 1)
example_gp_partab <- within(example_gp_partab, fixed[names == 't_switch'] <- 1)
example_gp_partab <- within(example_gp_partab, fixed[names == 'level_switch'] <- 1)
example_gp_partab <- within(example_gp_partab, fixed[names == 'prob_detect'] <- 1)


#################### MCMC setup #########################

#MCMC chain options
mcmc_pars <- c("iterations"=500000,"popt"=0.44,"opt_freq"=2000,
               "thin"=2500,"adaptive_period"=200000,"save_block"=100)
## Set pointer to the SEIR model as the incidence function
incidence_function <- gaussian_process_model

## Set standard deviations of prior distribution
sds_gp <- c("obs_sd"=0.5,"viral_peak"=2,
            "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
            "prob_detect"=0.03,
            "incubation"=0.25, "infectious"=0.5)

## Define a function that returns the log prior probability for a given vector of parameter
## values in `pars`, given the prior means and standard deviations set above.
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

for(horizon in horizons){
  ############### DATA #################
  ##############*******#################
  
  # Load in the Dental Data
  data <- read_csv(paste0(HOME_WD, data_dir, '/phase2_1_first_pos_msp.csv'), col_types=cols(Ct.e = "d", Ct.Rdrp = "d"))
  
  dataRAW <- read_csv(paste0(HOME_WD, data_dir, '/phase2_1_first_pos_msp.csv'), col_types=cols(Ct.e = "d", Ct.Rdrp = "d"))
  
  data <- data %>%
    mutate(collection_date = as.Date(collection_date, "%m/%d/%Y"))
  
  dataRAW <- dataRAW %>%
    mutate(collection_date = as.Date(collection_date, "%m/%d/%Y"))
  
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
           ct = Ct.e) %>%
    mutate(t = t+t_shift)%>%
    filter(!is.na(t))
  
  data <- data %>%
    mutate(t = calendar_time,
           ct = Ct.e) %>%
    filter(!is.na(t)) %>%
    select(t, ct)
  
  posdataRAW <- dataRAW %>%
    filter(ct != 40)
  
  # Uncomment me to use Negative Values
  data <- data%>%
    filter(ct != 40)
  
  
  data <- data %>%
    mutate(t = t+t_shift)
  
  
  data <- data %>%
    filter(t %in% horizon)
  
  
  ## This is for the GP version
  t_start <- beginning_times - 1
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
  
  posterior_function <- create_posterior_func(parTab=par_tab, 
                                              data=data, 
                                              PRIOR_FUNC=prior_func_gp, 
                                              INCIDENCE_FUNC=incidence_function,
                                              t_dist=t_dist, use_pos = TRUE)
  posterior_function(par_tab$values)
  
  
  dir.create("mcmc_chains/readme_multiple_cross_section",recursive=TRUE)
  ##################################
  ## RUN THE MCMC FRAMEWORK
  ## Run 3 MCMC chains. Note that it is possible to parallelize this loop with foreach and doPar
  ## Note the `use_pos` argument needs to be set here too
  # choose your horizons list here
  
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
                       filename=paste0("mcmc_chains/readme_multiple_cross_section/readme_gp_",toString(horizon), "_", chain_no),
                       CREATE_POSTERIOR_FUNC=create_posterior_func,
                       mvrPars=NULL,
                       OPT_TUNING=0.2,
                       use_pos=TRUE,
                       t_dist=t_dist)
  }
}
## Read in the MCMC chains
chains <- load_mcmc_chains(location="mcmc_chains/readme_multiple_cross_section",
                           parTab=par_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=TRUE,
                           unfixed=TRUE,
                           multi=FALSE)

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

# Create the horizons as an unlisted character to add to chains melted

# Add horizons to chains melted
chains_melted <- chains_melted %>%
  mutate(horizon_t = case_when(chain <= 3 ~ "horizon 1",
                               chain > 3 &  chain <= 6 ~ "horizon 2",
                               chain > 6 & chain <= 9  ~ "horizon 3"))
## Violin Plots for the horizons
violin <- ggplot(chains_melted) +
  geom_violin(aes(x=factor(horizon_t),y=value, fill=factor(horizon_t)), draw_quantiles = c(0.5)) +
  facet_wrap(~name, scales = "free_y")
violin

## Load in MCMC chains again, but this time read in the fixed parameters too 
## to ensure that the posterior draws are compatible with the model functions
chains <- load_mcmc_chains(location="mcmc_chains/readme_multiple_cross_section/",
                           parTab=par_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=TRUE,
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

x <- names(chain_comb)
is_dup <- stringr::str_detect(x,"prob")
names(chain_comb) <- if_else(is_dup,stringr::str_glue("{x}_{cumsum(is_dup)}"),x)

chain_comb <- chain_comb %>%
  mutate(horizon_t = case_when(chain <= 3 ~ "horizon 1",
                               chain > 3 &  chain <= 6 ~ "horizon 2",
                               chain > 6 & chain <= 9  ~ "horizon 3"))
## Load in true incidence curve to compare to our prediction
data(example_seir_incidence)
horizon_lst <- c('horizon 1', 'horizon 2', 'horizon 3')

names(chain_comb) <- if_else(stringr::str_detect(names(chain_comb),"prob"),"prob",names(chain_comb))

for(horizon in horizon_lst){
  chain_comb_filtered <- chain_comb %>%
    filter(horizon_t == horizon)
  predictions <- plot_prob_infection(chain_comb_filtered,nsamps=10, INCIDENCE_FUNC=incidence_function,
                                     solve_times=0:max(data$t),obs_dat=data,
                                     true_prob_infection=example_seir_incidence,
                                     smooth=TRUE) ## Smooth the trajectories a bit
  p_incidence_prediction <- predictions$plot + scale_x_continuous(limits=c(0,40))
  print(p_incidence_prediction)
}


# plot incidence vs prediction (NOT PROBABILITY)

# Load in Vital-E incidence curve
# reconstruct incidence from linelist
get_incidence_from_linelist <- function(data){
  data %>%
    group_by(t) %>%
    summarise(incidence = n())
}
posdataRAW <- dataRAW %>%
  filter(ct != 40)
incidence <- get_incidence_from_linelist(posdataRAW)

incidence <- incidence %>%
  mutate(prob_infection = incidence,
         t = t)

counter = 1
x <- names(chain_comb)
x[2] <- "overall_prob"
x[16] <- "prob_detect"
names(chain_comb) <- x
start_date <- min(posdataRAW$collection_date)

for(horizon in horizon_lst){
  x <- names(chain_comb)
  is_dup <- stringr::str_detect(x,"^prob$")
  names(chain_comb) <- if_else(is_dup,stringr::str_glue("{x}_{cumsum(is_dup)}"),x)
  chain_comb_filtered <- chain_comb %>%
    filter(horizon_t == horizon)
  
  
  
  names(chain_comb) <- if_else(stringr::str_detect(names(chain_comb),"prob_[:digit:]"),"prob",names(chain_comb))
  names(chain_comb_filtered) <- if_else(stringr::str_detect(names(chain_comb_filtered),"prob_[:digit:]"),"prob",names(chain_comb_filtered))
  horizon_time <- horizons[[counter]]
  obs_data <- filter_horizon_data(dataRAW, 0, horizon_time)
  predictions <- plot_infection_incidence(chain_comb_filtered, 
                                          nsamps=100,
                                          rawdata = dataRAW,
                                          posdata = posdataRAW,
                                          INCIDENCE_FUNC=incidence_function,
                                          solve_times=0:max(data$t),
                                          tshift = t_shift,
                                          obs_dat=obs_data,
                                          true_prob_infection=incidence, population = nrow(posdataRAW), horizons = FALSE, mode='Incidence',
                                          start_date = start_date)
  
  p_incidence_prediction <- predictions$plot 
  print(p_incidence_prediction)
  print("Printing R0 through chains")
  print(horizon_time)
  print(paste0("Median: ", median(chain_comb_filtered$R0)))
  print(paste0("Range: ", quantile(chain_comb_filtered$R0, .975) - quantile(chain_comb_filtered$R0, .025)))
  print(paste0("97.5%: " , quantile(chain_comb_filtered$R0, .975)))
  print(paste0("2.5%: " , quantile(chain_comb_filtered$R0, .025)))
  counter = counter + 1
}


# Repartition the horizons
selected_data <- data %>%
  filter(t %in% c(9,10,14,15,16))
## Use create_posterior_func to return the predicted Ct distribution rather than the posterior probability
model_func_gp <- create_posterior_func(par_tab,selected_data,NULL,incidence_function,"model")
## Pass model_func to a plotting function to observe predicted Ct distribution against data
p_distribution_fit_gp <- plot_distribution_fits(chain_comb, selected_data, model_func_gp,100,pos_only=TRUE)
## Joining, by = "t"
## Joining, by = c("t", "sampno")
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
p_distribution_fit_gp



