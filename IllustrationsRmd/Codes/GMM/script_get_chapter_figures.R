rm(list = ls())


# Graph parameters --------------------------------------------------------

graph_width <- 5
graph_height <- 5
graph_unit <- "in"

# Librairies --------------------------------------------------------------

library(tidyverse) # For data manipulation
library(GGally) # For ggpairs function
library(xtable) # To print teX tables
library(FactoMineR) # For PCA
library(factoextra)
library(traitor)


# Load data ---------------------------------------------------------------

# phytoplankton <- read.table("Phytoplankton.csv",
#                             sep = ",", header = TRUE, 
#                             row.names = 1)
donnees <- read.table("Bohemia_vegetation.csv",
                      header = TRUE, sep = ",") %>% 
  rename("Seed weight" = "Seed.weight")

# Data header ----------------------------------------------------------

xtable(head(donnees))

# Descriptive figures -----------------------------------------------------

## Pair plot ---------------------------------------------------------------

pair_plot <- ggpairs(donnees %>% select_if(is.numeric), 
                           upper = list(continuous = wrap("points", 
                                                          size = 1/.pt)),
        lower = "blank") +
  theme_bw() +
  theme(text = element_text(size = 10),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 
# ggsave("../../Figures/bohemia-pair-plot.png", 
#        plot = pair_plot, 
#        width = graph_width, 
#        height = graph_height, 
#        units = graph_unit)

## Load GMM results --------------------------------------------------------

if(!file.exists("gmm_est.RData")){
  source("script_get_gmm_estimates.R")
}
load("gmm_est.RData")

## Plot model selection graph ----------------------------------------------

criterions_plot <- map_dfr(gmm_est_multiple_ks, "final_criterions") %>% 
  mutate(log_lik = -2 * log_lik) %>% 
  rename("Negative log lik." = "log_lik") %>% 
  pivot_longer(-k) %>% 
  ggplot() + 
  aes(x = k, y = value, color = name) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = 1:6) +
  theme_bw() +
  labs(x = "Number of components",
       y = "Criterion value",
       color =  "")
# ggsave("../../Figures/bohemia-criterion-plot.png",
#        plot = criterions_plot,
#        width = graph_width,
#        height = graph_height,
#        units = graph_unit)

## Plot likelihoods evolution ----------------------------------------------

likelihoods_plot <- map(gmm_est_k3, "log_likelihoods") %>% 
  map_dfr(function(x){
    if(is.null(x)){
      return(NULL)
    }
    else{
      n_iter <- length(x) - 1
      return(data.frame(iteration = 0:n_iter,
                        log_lik = 2 * x * nrow(donnees)))
    }
  }, .id = "Replicate") %>% 
  ggplot() + 
  aes(x = iteration, y = log_lik, group = Replicate) +
  geom_line(alpha = .5) + 
  theme_bw() +
  lims(y = c(-2000, -1500)) +
  labs(x = "EM iteration", y = "Log likelihood value")

# ggsave("../../Figures/bohemia-likelihoods-plot.png",
#        plot = likelihoods_plot,
#        width = graph_width,
#        height = graph_height,
#        units = graph_unit)

## Plot clusters estimation ----------------------------------------------


estimated_clusters <- gmm_est_best$posterior_weights %>% 
  apply(1, which.max) %>% 
  factor()
cluster_pair_plot <- ggpairs(donnees %>% 
                       select_if(is.numeric) %>% 
                       mutate(Cluster = estimated_clusters), 
                     mapping = aes(color = Cluster, alpha = 1),
                     upper = list(continuous = wrap("points", alpha = 1,
                                                    size = 1/.pt)),
                     lower = "blank") +
  theme_bw() +
  theme(text = element_text(size = 10),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 

ggsave("../../Figures/bohemia-cluster-pair-plot.png",
       plot = cluster_pair_plot,
       width = graph_width,
       height = graph_height,
       units = graph_unit)

