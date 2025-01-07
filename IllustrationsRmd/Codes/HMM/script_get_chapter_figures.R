rm(list = ls())

library(tidyverse)
library(sf)
graph_width <- 5
graph_height <- 5
graph_unit <- "in"

# zone_map --------------------------------------------------------------

zone_map <- bind_rows(gadm("Canada", path = "./", level = 2) %>% 
                        st_as_sf() %>% 
                        filter(str_detect(NAME_1, "Newfoundland"), 
                               as.numeric(CC_2) %in% (1:9)) %>% # Only newfoundland without labrador
                        dplyr::select(COUNTRY), 
                      gadm("Saint Pierre and Miquelon", 
                           path = "./", level = 0) %>% 
                        st_as_sf() %>% 
                        dplyr::select(COUNTRY))

if(!file.exists("gannet_trajectory.csv")){
  source("script_create_clean_dataset.R")
}
donnees_brutes <- read.table("gannet_trajectory.csv", 
                      header = TRUE, sep = ",",
                      colClasses = c(t = "POSIXct")) 

library(xtable)  
donnees_brutes %>% 
  mutate(t = as.character(t)) %>% 
  head() %>% 
  rename("log10(Step length)" = log10_step) %>% 
  xtable(digits = 3) %>% 
  print(include.rownames = FALSE)



donnees <-   st_as_sf(donnees_brutes,
                      coords = c("Longitude", "Latitude"), crs = st_crs(zone_map)) 

trajectory_plot <- ggplot(zone_map) + 
  geom_sf(fill = "darkgreen") +
  theme_bw() +
  geom_sf(data = donnees %>% 
            st_combine() %>% 
            st_cast("LINESTRING")) +
  geom_sf(data = donnees, mapping = aes(color = t)) + 
  scale_color_viridis_c() + 
  theme(legend.position = "none") +
  coord_sf(xlim = c(-57, -52),
           ylim = c(46, 48))

ggsave("../../Figures/gannet-trajectory-plot.png",
       plot = trajectory_plot,
       width = graph_width,
       height = graph_height,
       units = graph_unit)

step_length_plot <- ggplot(donnees) + 
  aes(x = log10_step) +
  geom_histogram(fill = "lightblue", color = "black", 
                 mapping = aes(y = after_stat(density)), bins = 100) +
  theme_bw() +
  labs(y = "Empirical density", x = expression(log[10]~"of step length"))

ggsave("../../Figures/gannet-steplength-hist.png",
       plot = step_length_plot,
       width = graph_width,
       height = graph_height,
       units = graph_unit)


all_post_densities <- read.table("results_gannets_post_densities.txt",
                                 sep = ";", header = TRUE)

post_densities_plot <- ggplot(all_post_densities) + 
  aes(x = log10_step) + 
  geom_histogram(data = map_dfr(unique(all_post_densities$n_states),
                                function(k){
                                  mutate(donnees, n_states = k)
                                }),
                 mapping = aes(y = after_stat(density)), bins = 100,
                 color = "black", fill = "lightgray") +
  geom_line(aes(y = density, color = factor(state, levels = 1:10))) +
  facet_wrap(~n_states) +
  labs(color = "State", y = "Empirical density", 
       x = expression(log[10]~"of step length")) +
  theme_bw() 

ggsave("../../Figures/gannet-post_state_dist.png",
       plot = post_densities_plot,
       width = graph_width,
       height = graph_height,
       units = graph_unit)

model_criterions <- read.table("results_gannets_model_selection_criterions.txt",
                               sep = ";", header = TRUE)

model_criterions_plot <- ggplot(model_criterions,
       aes(x = K, y = value, color = criterion)) + 
  geom_point() + 
  geom_line() + 
  labs(x = "Number of hidden states", y = "Criterion value",
       color = "") +
  theme_bw()

ggsave("../../Figures/gannet-model-criterions.png",
       plot = model_criterions_plot,
       width = graph_width,
       height = graph_height,
       units = graph_unit)

illustrative_model <- readRDS("results_gannets_fitted_model.rds")

track_lines <- donnees %>% 
  mutate(
    viter = illustrative_model$best_seq,
    geometry_lagged = lead(geometry)
  ) %>%
  slice(-n()) %>% 
  # drop the NA row created by lagging
  mutate(
    line = st_sfc(purrr::map2(
      .x = geometry, 
      .y = geometry_lagged, 
      .f = ~{st_union(c(.x, .y)) %>% st_cast("LINESTRING")}
    ), crs = st_crs(.)))

trajectory_plot_with_states <- ggplot(zone_map) + 
  geom_sf(fill = "darkgreen") +
  theme_bw() +
  geom_sf(data = track_lines, aes(geometry = line, color = factor(viter))) + 
  geom_sf(data = donnees %>% 
            mutate(viter = illustrative_model$best_seq %>% factor()), 
          mapping = aes(color = viter)) + 
  theme(legend.position = "inside", 
        legend.direction = "horizontal", legend.background = element_blank(),
        legend.position.inside = c(0.4, 0.9)) +
  coord_sf(xlim = c(-56.3, -56),
           ylim = c(46.75, 46.95)) + 
  # scale_color_manual(values = ten_color_palette[1:6]) + 
  labs(color = "Speed regime")

ggsave("../../Figures/gannet-trajectory-plot-with-states.png",
       plot = trajectory_plot_with_states,
       width = graph_width,
       height = graph_height,
       units = graph_unit)

transition_matrix_plot <- ggcorrplot(orig_mat[, 6:1], 
                                     lab = TRUE, 
                                     type = "full", lab_col = "red") +
  theme(legend.position = "none") +
  scale_fill_viridis_c()
ggsave("../../Figures/gannet-transition-mat.png",
       plot = transition_matrix_plot,
       width = graph_width,
       height = graph_height,
       units = graph_unit)
orig_mat <- illustrative_model$P
colnames(orig_mat) <- rownames(orig_mat) <- paste0("S", 1:ncol(orig_mat))
colnames(my_mat)
rownames(my_mat) <- nrow(illustrative_model$P):1


corrplot::corrplot(illustrative_model$P, col.lim = c(0, 1))
ggcorrplot(illustrative_model$P) +
  scale_fill_gradientn(colors = viridis(256, option = 'D'))
