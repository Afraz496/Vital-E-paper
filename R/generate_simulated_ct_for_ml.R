# Generate Simulations for Ct
library(virosolver)
library(tidyverse)
library(patchwork)
library(lazymcmc)
library(foreach)
library(doParallel)
library(here)

args <- commandArgs(trailingOnly = TRUE)

print(args)
print(typeof(args[2]))
# CHANGE TO 10, 100, 1000
N <- 500000
print(N)
sample_days <- seq(1,140,by=1)
sample_size <- args[2]
set.seed(args[1])


data(example_seir_partab)
pars <- example_seir_partab$values
names(pars) <- example_seir_partab$names
pars["t0"] <- 0

times <- seq(1,200,by=1)

## Pull parameters for SEIR model
seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["incubation"],1/pars["infectious"])
## Set up initial conditions.
## Note if population_n=1, then solves everything per capita
init <- c((1-pars["I0"])*N,0,pars["I0"]*N,0,0,0)

## Solve the SEIR model using the rlsoda package
#sol <- rlsoda::rlsoda(init, times, C_SEIR_model_rlsoda, parms=seir_pars, dllname="virosolver",
#                      deSolve_compatible = TRUE,return_time=TRUE,return_initial=TRUE,atol=1e-10,rtol=1e-10)

## Solve the SEIR model using the lsoda package. lsoda runs about 4x slower than rlsoda, but
## the lsoda package is available on CRAN, making it more user-friendly.
sol <- deSolve::ode(init, times, func="SEIR_model_lsoda",parms=seir_pars,
                    dllname="virosolver",initfunc="initmodSEIR",
                    nout=0, rtol=1e-6,atol=1e-6)

## Convert to data frame and column names
sol <- as.data.frame(sol)
colnames(sol) <- c("time","S","E","I","R","cumu_exposed","cumu_infected")
## Get Rt
## 
sol$Rt <- (sol$S) * pars["R0"] 

## Shift time for start
sol$time <- sol$time + floor(pars["t0"])

## Dummy rows from pre-seeding
if(pars["t0"] > 0){
  dummy_row <- data.frame("time"=0:(floor(unname(pars["t0"]))-1),"S"=N,"E"=0,"I"=0,"R"=0,"cumu_exposed"=0,"cumu_infected"=0,"Rt"=unname(pars["R0"])*N)
  sol <- bind_rows(dummy_row, sol)
}
sol <- sol[sol$time %in% times,]

## Pull raw incidence in absolute numbers
inc <- c(0,diff(sol[,"cumu_exposed"]))
inc <- pmax(0, inc)

## Get per capita incidence and prevalence
per_cap_inc <- inc/N
per_cap_prev <- (sol$E + sol$I)/N

## Melt solution (wide to long format) and get per capita
sol <- reshape2::melt(sol, id.vars="time")
sol$value <- sol$value/N






prob_infection <- sol %>% filter(variable == "I") %>% pull(value)
infection_times <- simulate_infection_times(N,prob_infection,overall_prob = 1)


# get sample day
sample_t <- 100
# sample 100 from population:
sample_from_population <- sample(infection_times,100)
times_since_infected <- if_else(sample_t > sample_from_population, 
                                sample_t - sample_from_population,
                                -1
)
tibble(y=times_since_infected) %>%
  ggplot(aes(x="Sample Time",y=y)) +
  geom_violin() +
  geom_jitter()

sampled_times_since_infected <- tibble()
for(sample_t in sample_days){
  # sample 100 from population:
  sample_from_population <- sample(infection_times,sample_size)
  times_since_infected <- if_else(sample_t > sample_from_population, 
                                  sample_t - sample_from_population,
                                  -1
  )
  sampled_times_since_infected <- sampled_times_since_infected %>% 
    bind_rows(
      tibble(times_since_infected = times_since_infected,
             sample_t = sample_t, Rt = sol %>% filter(variable == "Rt",time == {{sample_t}}) %>% pull(value))
    )
}



# viral_loads <- simulate_viral_loads_example(sampled_times_since_infected$times_since_infected,
#                                             pars,N=1)

times_since_infected <- sampled_times_since_infected$times_since_infected

t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
sd_mod <- rep(pars["sd_mod"], max(times_since_infected))
unmod_vec <- 1:min(t_switch,max(times_since_infected))
sd_mod[unmod_vec] <- 1
decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
## For the next sd_mod_wane days, variance about Ct trajectory decreases linearly
sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])

sampled_ct <- matrix(ncol=1,nrow=length(times_since_infected))
for(i in 1:length(times_since_infected)){
  age <- times_since_infected[i]
  ## For each value we're going to simulate, pre-determine if it will still be detectable or not
  detectable_statuses <- rnbinom(1, 1, prob=pars["prob_detect"]) + 
    pars["tshift"] + pars["desired_mode"] + pars["t_switch"]
  cts <- rep(virosolver::viral_load_func(pars, age, FALSE, 0),1)
  cts[detectable_statuses <= age] <- 1000
  sd_used <- pars["obs_sd"]*sd_mod[age]
  ## Generate a observations of Ct values from gumbel distribution for a specified mode
  ct_obs_sim <- extraDistr::rgumbel(1, cts, sd_used)
  ## Set Ct values greater than intercept to intercept value
  ct_obs_sim <- pmin(ct_obs_sim, pars["intercept"])
  sampled_ct[i,] <- ct_obs_sim
}
sampled_ct <- data.frame(sampled_ct)

sampled_ct$time_since_infection <- times_since_infected



sampled_ct <- sampled_times_since_infected %>%
  bind_cols(sampled_ct %>% select(-time_since_infection)) 




sampled_ct %>% 
  rename("time" = "sample_t","ct" = "sampled_ct") %>%
  select(time,ct,Rt) %>% 
  write_csv(here(paste0("Simulation Files Sizes/simulated_ct_output_",args[1],"_",args[2],".csv")))

