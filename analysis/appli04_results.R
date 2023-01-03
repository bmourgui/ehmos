############################################%
# Results analysis presented in the manuscript
#
# last modification: 16/12/21 (create file)
# bastien.mourguiart@gmail.com
#
# file name: appli04_results.R
#
############################################%


#### Averaging results ----
appli_res  %>%
  dplyr::group_by(model) %>%
  dplyr::summarise_at("shift", .funs = c(mean, sd)) %>%
  dplyr::rename("m.shift"="fn1", "sd.shift"="fn2") -> appli_mean_shift

appli_res %>%
  dplyr::mutate("CI.shift" = upr.shift - lwr.shift) %>%
  dplyr::group_by(model) %>%
  dplyr::summarise_at("CI.shift", .funs = c(mean, sd)) %>%
  dplyr::rename("m.CI.shift"="fn1", "sd.CI.shift"="fn2") -> appli_mean_CI.shift


#### Figure 3: Plot species shift estimates ----

## create shortcuts for species names
appli_res %>%
  tidyr::separate(sp, c("genre", "species")) %>%
  dplyr::mutate("genre" = substr(genre, 1, 2),
                "sp" = paste(genre, species, sep = ". ")) -> res_graph

## Diminished credible intervals for visualization because GLMM intervals could be too huge and render difficult graphic readability
for (i in 1:nrow(res_graph)){
  res_graph$lwr.shift[i] <- max(res_graph$lwr.shift[i], -800)
  res_graph$upr.shift[i] <- min(res_graph$upr.shift[i], 800)
}


# Make a vector with ordered species according to their optimum positions
data.frame("species" = unique(res_graph$sp),
           "opti" = res_graph[res_graph$model == "ehmos",]$opt) %>%
  dplyr::arrange(opti) %>%
  dplyr::select(species) -> species_ordered

species2 <- factor(res_graph$sp,
                   level = as.character(species_ordered$species))


## Look at species having a mean estimated width larger than 1200 m (which is the minimal value in simulated generalist species)
appli_res %>%
  dplyr::mutate("width" = (width1 + width2) / 2) %>%
  dplyr::select(sp, model, width) %>%
  tidyr::spread(model, width) %>%
  dplyr::mutate("m.width" = (ehmos + glmm + ttest)/3,
                "is.generalist" = as.factor(m.width > 1200)) -> is_generalist_d

appli_res %>%
  dplyr::mutate("width" = (width1 + width2) / 2) %>%
  dplyr::select(sp, model, width) %>%
  dplyr::mutate("sp" = factor(res_graph$sp,
                              level = as.character(species_ordered$species))) %>%
  tidyr::spread(model, width) %>%
  dplyr::mutate("m.width" = (ehmos + glmm + ttest)/3) %>%
  dplyr::filter(m.width > 1200) %>%
  dplyr::select(sp) -> list.generalist

levels(is_generalist_d$is.generalist) <- c("black","purple")

for (i in 1:length(levels(species2))){
  if (levels(species2)[i] %in% list.generalist$sp){
    levels(species2)[i] <- paste0(levels(species2)[i], "*")
  }
}

## Create coloured rectangles to distinguish middle from edge species
rects_app <- data.frame(xstart = c(0,9.5, 22.5),
                        xend = c(9.5,22.5, 24.5),
                        col = c("orange","forestgreen","orange"))


# Figure: Species specific shifts estimates by model
res_graph %>%
  ggplot() +
  geom_errorbar(aes(x = species2,
                    ymin = lwr.shift,
                    ymax = upr.shift,
                    colour = model),
                position = position_dodge(width = 0.8)
                ) +
  geom_point(aes(x = species2,
                 y = shift,
                 colour = model,
                 shape = model),
             position = position_dodge(width = 0.8),
             size = 2
             ) +
  geom_vline(xintercept = seq(1.5, 23.5, 1), 
             color = "grey",
             lty = 3,
             size = 0.3) +
  ylim(c(-800, 800)) +
  theme_bw() +
  theme(
    #legend.position = c("top"),
    legend.box.background = element_rect(color="black"),
    text = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(
          size = 8,
          angle = -60,
          hjust = 0,
          colour = rep(c("orange", "forestgreen", "orange"), c(9,13,2))
        ),
        axis.title.x = element_blank(),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
        axis.ticks.x = element_line(colour = "grey")) +
  ylab("Estimated shift") +
  geom_hline(yintercept = c(-800, 800), lty = 1)  +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02"),
                     name = "Model",
                     labels = c("EHMOS", "cGLMM", "t-test")) +
  scale_shape_manual(values = c(19, 17, 15),
                     name = "Model",
                     labels = c("EHMOS", "cGLMM", "t-test")) +
  scale_fill_manual(values = c("forestgreen", "orange"),
                    name = "Optimum type",
                    labels = c("Middle", "Edge")) -> fig3

ggsave(filename = here::here("results", "figs", "fig3.png"),
       plot = fig3,
       height = 140,
       width = 168,
      units = "mm")

#### Paired t-test between method estimates ----
t.test(appli_res[appli_res$model == "ttest", ]$shift,
       appli_res[appli_res$model == "glmm", ]$shift,
       paired = TRUE) -> appli_ttest.vs.glmm

t.test(appli_res[appli_res$model == "ttest", ]$shift,
       appli_res[appli_res$model == "ehmos", ]$shift,
       paired = TRUE) -> appli_ttest.vs.ehmos


#### Goodness-of-fit (Bayesian p-value) - EHMOS ----
run.gof <- FALSE
if (run.gof == TRUE){
  plot(x = out_appli_ehmos$sims.list$p.fit, 
       y = out_appli_ehmos$sims.list$p.fitnew,
       xlab = "p.fit",
       ylab = "p.fitnew",
       xlim = range(c(out_appli_ehmos$sims.list$p.fit, out_appli_ehmos$sims.list$p.fitnew)),
       ylim = range(c(out_appli_ehmos$sims.list$p.fit, out_appli_ehmos$sims.list$p.fitnew)),
       main = "All species")
  abline(a = 0, b = 1)
  par(mfrow = c(4, 6))
  for (i in 1:24){
    pfit <- out_appli_ehmos$sims.list$p.fitSP[, i]
    pfitnew <- out_appli_ehmos$sims.list$p.fitnewSP[, i]
    bp <- sum(pfit > pfitnew)/length(pfit)
    plot(x = pfit, 
         y = pfitnew,
         xlab = "p.fit",
         ylab = "p.fitnew",
         xlim = range(c(pfit, pfitnew)),
         ylim = range(c(pfit, pfitnew)),
         main = paste0("species ", i))
    mtext(text = paste0("bp = ", round(bp, 2)),
          cex = 0.8)
    abline(a = 0, b = 1)
  }
  
  
  #### Goodness-of-fit (Bayesian p-value) - GLMM ----
  plot(x = out_appli_glmm$sims.list$p.fit, 
       y = out_appli_glmm$sims.list$p.fitnew,
       xlab = "p.fit",
       ylab = "p.fitnew",
       xlim = range(c(out_appli_glmm$sims.list$p.fit, out_appli_glmm$sims.list$p.fitnew)),
       ylim = range(c(out_appli_glmm$sims.list$p.fit, out_appli_glmm$sims.list$p.fitnew)),
       main = "All species")
  abline(a = 0, b = 1)
  par(mfrow = c(4, 6))
  for (i in 1:24){
    pfit <- out_appli_glmm$sims.list$p.fitSP[, i]
    pfitnew <- out_appli_glmm$sims.list$p.fitnewSP[, i]
    bp <- sum(pfit > pfitnew)/length(pfit)
    plot(x = pfit, 
         y = pfitnew,
         xlab = "p.fit",
         ylab = "p.fitnew",
         xlim = range(c(pfit, pfitnew)),
         ylim = range(c(pfit, pfitnew)),
         main = paste0("species ", i))
    mtext(text = paste0("bp = ", round(bp, 2)),
          cex = 0.8)
    abline(a = 0, b = 1)
  }
}

