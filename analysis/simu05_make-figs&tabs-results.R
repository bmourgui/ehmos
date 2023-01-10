


#### 1. Figure 1: Simulation workflow ----
## .. 1.1. Define graphic elements that are used multiple times ####
theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(size = 8),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"))) -> th # general theme

rects <- data.frame(xstart = c(1000,
                               quantile(X_A1, 0.1),
                               quantile(X_A1, 0.9),
                               1000,
                               quantile(X_A2, 0.1),
                               quantile(X_A2, 0.9)
                               ),
                    xend = c(quantile(X_A1, 0.1),
                             quantile(X_A1, 0.9),
                             3000,
                             quantile(X_A2, 0.1),
                             quantile(X_A2, 0.9),
                             3000),
                    col = c("orange","grey","orange"),
                    scenario = rep(c("A2","A1"), each = 3))

geom_rect(data = rects,
          aes(xmin = xstart,
              xmax = xend,
              ymin = -Inf,
              ymax = Inf,
              fill = col),
          alpha = 0.2)  -> grect # add colored rectangles for distinguish the middle and edges of sampling range

rcols <- c("forestgreen", "orange") # colors for rectangles
#####                                                 ####%

## .. 1.2. Site repartition ####
alti_site <- data.frame(scenario = rep(c("A1","A2"), each = 300),
                        elevation = c(X_A2, X_A1))

ggplot(alti_site) +
  geom_histogram(aes(x = elevation),
                 bins = 10,
                 col = "black",
                 fill = "grey80") +
  xlab("Elevation (m)") +
  ylab("Number of sites") +
  facet_grid(. ~ scenario) +
  xlim(c(1000,3000)) +
  scale_color_manual(values = rcols) +
  scale_fill_manual(values = rcols) +
  th +
  grect -> g_subscenarA


## .. 1.3. Optima repartition ####
sim.data$scenario <- as.factor(sim.data$scenario)
levels(sim.data$scenario) <- c("A2xB1xC1","A2xB1xC2","A2xB2xC1","A2xB2xC2","A1xB1xC1","A1xB1xC2","A1xB2xC1","A1xB2xC2")

sim.data %>%
  dplyr::filter(scenario %in% c("A1xB1xC1", "A2xB2xC2")) %>%
  dplyr::mutate("scenario" = substr(scenario, 1,2)) %>%
  ggplot() +
  geom_vline(xintercept = rep(1:20, 2), lty=2, color = "grey90") + # manual panel grid
  geom_hline(yintercept = seq(800,3200,200), lty=2, color = "grey90")  + # manual panel grid
  geom_point(aes(x = rep(1:20, 2), y = Opti, color = T_opti)) + # Optima at first sampling occasion
  geom_point(aes(x = rep(1:20, 2), y = (Opti + Shift), color = T_opti), shape = 1) + # Optima after shift
  coord_flip() +
  xlab("Species") +
  ylab("Elevation (m)") +
  facet_grid(~ scenario,
             labeller = labeller(scenario = c("A1" = "B1", "A2" = "B2"))
             ) +
  scale_color_manual(values = rcols) +
  scale_fill_manual(values = rcols) +
  th +
  geom_rect(data = rects,
            aes(ymin = xstart,
                ymax = xend,
                xmin = -Inf,
                xmax = Inf,
                fill = col),
            alpha = 0.2) -> g_subscenarB


## .. 1.4. Species response curves (specialist vs generalist) ####
## p and w used for simulate data
psi_S <- SRC(pm = 0.925,
             pinf = .85,
             psup = 0.99,
             wm = 750,
             winf = 500,
             wsup = 900,
             opt = 2000,
             x = 1000:3000)

psi_G <- SRC(pm = 0.725,
             pinf = .65,
             psup = 0.8,
             wm = 1450,
             winf = 1200,
             wsup = 1600,
             opt = 2000,
             x = 1000:3000)

#plot the "range" of species response curves
data.frame("altitude" = 1000:3000,
           "m.psi" = c(psi_S[, 1], psi_S[, 1], psi_G[, 1]),
           "inf.psi" = c(psi_S[, 2], psi_S[, 2], psi_G[, 2]),
           "sup.psi" = c(psi_S[, 3], psi_S[, 3], psi_G[, 3]),
           "scenario" = rep(c("C1", "C2", "C2"), each = 2001),
           "sp" = rep(c("spe", "spe", "gen"), each = 2001)) %>%
  ggplot(aes(x = altitude, color = sp)) +
  geom_line(aes(y = m.psi)) +
  geom_line(aes(y = inf.psi), lty = 2) +
  geom_line(aes(y = sup.psi), lty = 2) +
  xlab("Elevation (m)") +
  ylab("Occupancy probability") +
  scale_color_manual(values = c("purple", "black")) +
  facet_grid(. ~ scenario) +
  th -> g_subscenarC


## .. 1.5. Species occupancy resulting of sub-scenario combination ####
reshape2::melt(psi,
               varnames = c("site", "sp", "period", "scenario")) %>%
  dplyr::mutate("alti" = c(rep(X_A1,  20*4*2), rep(X_A2,  20*4*2))) %>%
  dplyr::arrange(sp, period, scenario, alti) -> occu

occu$scenario <- as.factor(occu$scenario)
levels(occu$scenario) <- c("A2xB1xC1","A2xB1xC2","A2xB2xC1","A2xB2xC2","A1xB1xC1","A1xB1xC2","A1xB2xC1","A1xB2xC2")

occu %>%
  dplyr::filter(scenario %in% c("A1xB1xC1", "A2xB2xC2")) %>%
  dplyr::inner_join(sim.data, by=c("scenario", "sp")) %>%
  dplyr::mutate("scenario" = substr(scenario, 1, 2)) %>%
  ggplot(aes(fill = sp)) +
  geom_hline(yintercept = seq(0, 1, 0.2), lty=2, color="grey90") +
  geom_hline(yintercept = seq(0, 1, 0.2), lty=2, color="grey90") +
  geom_line(aes(x = alti, y = value, color = T_shape, linetype = period)) +
  scale_fill_manual(values = c(rcols, rep("black",20))) +
  scale_color_manual(values=c("purple", "black")) +
  facet_grid(. ~ scenario,
             labeller = labeller(scenario = c("A1" = "A1xB1xC1",
                                              "A2" = "A2xB2xC2")
                                 )
             ) +
  xlab("Elevation (m)") +
  ylab("Occupancy probability") +
  th +
  grect -> g_scenar

## .. 1.6. Final graph ####
cowplot::plot_grid(g_subscenarA, g_subscenarB, g_subscenarC, g_scenar, 
                   nrow=4,
                   labels=c("(a)", "(b)","(c)","(d)"),
                   label_size = 8) -> fig1

## .. 1.7. Save figure 1 ####
ggsave(filename = here::here("results", "figs", "fig1.png"),
       plot = fig1,
       height = 190,
       width = 168,
       units = "mm")




#### 2. Figure 2: Performance metrics by scenario ----
## .. 2.1. RMSE ####
perf_sp %>% # created in simu04_compute-performance-metrics.R
  ggplot(aes(x = model, y = RMSE)) %>%
  ggnested(ylab = "RMSE (m)") %>% # home made functions ggnested and ggtheme
  ggtheme() -> g.RMSE

## .. 2.2. Bias ####
perf_sp %>%
  ggplot(aes(x = model, y = bias)) %>%
  ggnested(ylab = "Bias (m)") %>%
  ggtheme() -> g.bias

## .. 2.3. Interval coverage ####
perf_sp %>%
  ggplot(aes(x = model, y = coverage)) %>%
  ggnested(ylab = "Interval coverage") %>%
  ggtheme() -> g.coverage

## .. 2.4. Interval score ####
perf_sp %>%
  ggplot(aes(x = model, y = IS)) %>%
  ggnested(ylab = "Interval score") %>%
  ggtheme() -> g.IS

## .. 2.5. Final graph ####
ggpubr::ggarrange(g.RMSE,
                  g.bias,
                  g.coverage,
                  nrow = 3) -> fig2

## .. 2.6. Save figure 2 ####
ggsave(filename = here::here("results", "figs", "fig2.png"),
       plot = fig2,
       height = 190,
       width = 168,
       units = "mm")


#### 3. Table 1: Interval scores by scenario ----

perf_sp %>%
  average_metrics(scenario, model) %>%
  dplyr::mutate("m.IS" = round(m.IS, 2),
                "sd.IS" = round(sd.IS, 2),
                "IS" = paste(m.IS, " (", sd.IS, ")",sep = "")) %>%
  dplyr::select(scenario, model, IS) %>%
  tidyr::spread(model, IS) %>%
  dplyr::bind_cols("pos" = c(5:8, 1:4)) %>% # manually change scenario order
  dplyr::arrange(pos) %>% # uniform sampling goes first
  dplyr::select(-pos) %>% # delete created column for changing order
  dplyr::rename("Scenario" = "scenario",
                "EHMOS" = "ehmos",
                "cGLMM" = "glmm",
                "t-test" = "t.test") -> tab1

save(tab1,
     file = here::here("results", "tables", "tab1.RData"))


#### 4. Table S1: Convergence issues by scenario ----

full.res.shift %>%
  dplyr::mutate("convergence_fail" = converge > 1) %>%
  dplyr::group_by(model, scenario, replication) %>%
  dplyr::summarise("convergence_fail" = sum(convergence_fail) > 1) %>%
  dplyr::group_by(model, scenario) %>%
  dplyr::summarise("convergence_fail" = sum(convergence_fail)) %>%
  tidyr::spread(model, convergence_fail) %>%
  dplyr::select(-t.test) %>%
  dplyr::rename("Scenario" = "scenario",
                "EHMOS" = "ehmos",
                "cGLMM" = "glmm") -> tabS1

save(tabS1,
     file = here::here("results", "tables", "tabS1.RData"))


## 5. Table S2: computation times ----
## Average by model (used in results)
full.res.shift %>%
  dplyr::group_by(model) %>%
  dplyr::summarise("mean" = round((mean(time) / 60), 0),
                   "SD" = round((sd(time) / 60), 0)) -> m.compute_time

## Table
full.res.shift %>%
  dplyr::group_by(scenario, model) %>%
  dplyr::summarise("time" = round((mean(time) / 60), 1)) %>%
  tidyr::spread(model, time) %>%
  dplyr::select(-t.test) %>%
  dplyr::rename("Scenario" = "scenario",
                "EHMOS" = "ehmos",
                "cGLMM" = "glmm") -> tabS2

save(tabS2,
     file = here::here("results", "tables", "tabS2.RData"))

## Table S3: All performance metrics by scenario ----

perf_sp %>%
  average_metrics(scenario, model) %>%
  dplyr::mutate_if(is.numeric, function(x)round(x, 2)) %>%
  dplyr::mutate(
    "RMSE" = paste(m.RMSE, " (", sd.RMSE, ")",sep = ""),
    "Bias" = paste(m.bias, " (", sd.bias, ")",sep = ""),
    "Interval coverage" = paste(m.coverage, " (", sd.coverage, ")",sep = ""),
    "IS" = paste(m.IS, " (", sd.IS, ")",sep = "")
    ) %>%
  dplyr::select(scenario, model, RMSE, Bias, `Interval coverage`, IS) %>%
  tidyr::gather("metric", "value", 3:6) %>%
  tidyr::spread(model, value) %>%
  dplyr::arrange(metric) %>%
  dplyr::select(- metric) %>%
  dplyr::rename("Scenario" = "scenario",
                "EHMOS" = "ehmos",
                "cGLMM" = "glmm",
                "t-test" = "t.test") -> tabS3


save(tabS3,
     file = here::here("results", "tables", "tabS3.RData"))


## Figure S2: Interval score ----
perf_sp %>%
  dplyr::filter(model != "glmm") %>%
  ggplot(aes(x = model, y = IS)) %>%
  ggnested(ylab = "Interval score") %>%
  ggtheme() -> g.IS2

figS2 <- ggpubr::ggarrange(g.IS, g.IS2, nrow = 2)

ggsave(filename = here::here("results", "figs", "figS2.png"),
       plot = figS2,
       height = 8.3,
       width = 11.7)
