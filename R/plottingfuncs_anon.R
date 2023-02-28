library(virosolver)
library(tidyverse)
library(patchwork)
library(lazymcmc)
library(foreach)
library(doParallel)
library(cowplot)
library(lubridate)
## Functions

# Return confidence interval to use for prevalence plots
ci_estimate_lower<- function(p1,n1){
  prop_est<- p1
  prop_std<-sqrt((p1*(1-p1))/n1)
  lci<- prop_est- 1.96*prop_std
  return(lci)
}

ci_estimate_higher <- function(p1,n1){
  prop_est<- p1
  prop_std<-sqrt((p1*(1-p1))/n1)
  uci <- prop_est + 1.96*prop_std
  return(uci)
}

plot_infection_incidence_anonymous <- function(chain,
                                               nsamps,
                                               posdata, rawdata,  
                                               INCIDENCE_FUNC,
                                               solve_times,
                                               obs_dat=NULL,
                                               true_prob_infection=NULL,
                                               tshift=0,
                                               smooth=FALSE, population = NULL, horizons = FALSE, mode='Incidence',
                                               start_date = as.Date("12-12-31", "%m-%d-%Y")){
  
  ## Take n samples from the MCMC chain
  samps <- sample(unique(chain$sampno),nsamps)
  all_res <- NULL
  for(i in seq_along(samps)){
    samp <- samps[i]
    ## Return parameters for each sample according to sample number
    tmp_pars <- lazymcmc::get_index_pars(chain, samp)
    
    ## Solve compartmental model (SEIR or SEEIRR)
    prob_infection_tmp <- INCIDENCE_FUNC(tmp_pars, solve_times)
    
    if(!is.null(population)){
      prob_infection_tmp <- prob_infection_tmp*population
    }
    
    ## Smooth out curve for plotting  
    if(smooth){
      prob_infection_tmp <- pmax(smooth.spline(prob_infection_tmp,spar=0.3)$y,0.0000001)
    }
    all_res[[i]] <- tibble(t=solve_times+tshift,prob_infection=prob_infection_tmp,sampno=i)
  }
  
  
  
  ## Combine results for all samples
  posterior_dat <- do.call("bind_rows",all_res)
  
  ## Return parameters with highest likelihood
  best_pars <- lazymcmc::get_best_pars(chain)
  
  ## Solve compartmental model using parameters with highest likelihood
  best_prob_infection <- INCIDENCE_FUNC(best_pars, solve_times)
  
  ## Smooth out curve for plotting  
  if(smooth){
    best_prob_infection <- pmax(smooth.spline(best_prob_infection)$y,0.0000001)
  }
  best_prob_dat <- tibble(t=solve_times+tshift,prob_infection=best_prob_infection,sampno="MAP")
  
  # Fix MAP for Incidence
  if(!is.null(population)){
    best_prob_dat$prob_infection <- best_prob_dat$prob_infection*population
  }
  get_incidence_from_linelist <- function(data){
    data %>%
      group_by(t) %>%
      summarise(incidence = n())
  }
  # Add the incidence of true data to figure
  posdata <- posdata %>%
    get_incidence_from_linelist()
  
  rawdata <- rawdata %>%
    get_incidence_from_linelist()
  
  combined_data <- rawdata %>%
    rename(total_sample = incidence) %>%
    inner_join(posdata,by="t") %>%
    mutate(prevalence = incidence/total_sample)
  
  
  
  # Add confidence intervals of the prevalence
  combined_data <- combined_data %>%
    mutate(low = ci_estimate_lower(prevalence,total_sample),
           high = ci_estimate_higher(prevalence,total_sample),
           collection_date = start_date + days(t))
  
  posterior_dat <- posterior_dat %>%
    mutate(collection_date = start_date + days(t))
  
  posdata <- posdata %>%
    mutate(collection_date = start_date + days(t))
  
  
  obs_dat <- obs_dat %>%
    mutate(collection_date = start_date + days(t))
  
  best_prob_dat <- best_prob_dat %>%
    mutate(collection_date = start_date + days(t))
  # Create the date time labels and breaks
  dates <- posdata %>%
    select(collection_date) %>%
    distinct() %>%
    mutate(labels = format(collection_date, "%b %d"))
  
  posterior_dat_grouped <- posterior_dat %>%
    group_by(collection_date,t) %>% 
    summarize(p50 = quantile(prob_infection, 0.5),
              p05 = quantile(prob_infection, 0.05),
              p95 = quantile(prob_infection, 0.95))
  
  breaks <- dates$collection_date
  
  labels <- dates$labels
  # Clean xticks in weekly format
  
  labels[seq(2, floor(0.25*length(labels)) - 1,1)] = ""
  labels[seq(floor(0.25*length(labels))+1, (floor(0.5*length(labels)) - 1),1)] = ""
  labels[seq(floor(0.5*length(labels))+1, (floor(0.75*length(labels)) - 1),1)] = ""
  labels[seq(floor(0.75*length(labels))+1, length(labels) - 1,1)] = ""
  
  ## Begin plotting
  if(mode == 'Prevalence'){
    p1 <- ggplot(posterior_dat) +
      ## Plot probability of infection for each posterior sample
      geom_line(aes(x=t,y=prob_infection,group=sampno),size=0.1, alpha = 0.3) +
      ## Plot one line for the MAP (maximum posterior probability of infection). This is the trajectory with the 
      ## highest likelihood.
      geom_line(data = combined_data, aes(x = t,y = prevalence), color = "Green") +
      geom_line(data=best_prob_dat,aes(x=t,y=prob_infection,col="Mean"),size=0.5)+
      geom_ribbon(data = combined_data, aes(ymin = low, ymax = high, x = t), fill = "grey70", alpha = 0.3) +
      geom_bar(data = posdata, aes(x = t), alpha= 0.3)+
      xlab("Calendar Time") +
      ylab("Prevalence") +
      theme_classic()
  }
  else if(mode == 'Incidence'){
    coeff = 1
    p1 <- ggplot(posterior_dat) +
      ## Plot probability of infection for each posterior sample
      geom_line(aes(x=t,y=prob_infection,group=sampno),size=0.1, alpha = 0.3) +
      geom_ribbon(data = posterior_dat_grouped, aes(ymin = p05, ymax = p95, x = t), fill = "Blue", alpha = 0.2) +
      geom_col(data = posdata, aes(x = t, y = incidence), fill= '#FF800E',alpha= 0.3)+
      scale_y_continuous(
        name = "Incidence",
        sec.axis = sec_axis(~. / coeff, name = "Number of Cases", )) +
      xlab("Time (in days)") +
      ylab("Incidence") +
      theme_classic()
  }
  ## Add the Median Chain as a green line to the plot with a CI
  
  
  ## Add vertical line with the sample date to the plot
  if(!is.null(obs_dat)){
    p1 <- p1 +
      geom_vline(data=data.frame(x=unique(obs_dat$t)),
                 aes(xintercept=x,linetype="Sample date"),
                 col="red",size=0.25)+
      geom_rect(data = data.frame(x=unique(obs_dat$t) ),
                aes(xmin = min(posdata$t), xmax = max(obs_dat$t), ymin = 0, ymax = Inf,
                    fill = "Observed Ct Values")
                , alpha = 0.1) 
  }
  
  ## Adding colors to the plot and modifying the legend
  p1 <- p1 +
    scale_color_manual(values=c("Posterior draw"="gray50")) +
    scale_linetype_manual(values=c("Sample date"="dashed")) +
    scale_fill_manual(NULL, 
                      values = 'gray70',
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    guides(col=guide_legend(title=NULL),linetype=guide_legend(title=NULL)) +
    theme(legend.position="bottom")
  return(list(predictions=posterior_dat, map_prediction=best_prob_dat, plot=p1))
}