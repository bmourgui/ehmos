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

levels(is_generalist_d$is.generalist) <- c("black","purple")

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
  geom_rect(data = rects_app,
            aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col),
            alpha = 0.4) +
  ylim(c(-800, 800)) +
  theme_bw() +
  theme(axis.text.x = element_text(
                        size = 16,
                        angle = -60,
                        hjust = 0,
                        colour = as.character(is_generalist_d$is.generalist)
                        ),
        axis.title.x = element_blank()) +
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
       height = 8.3,
       width = 11.7)

#### Paired t-test between method estimates ----
t.test(appli_res[appli_res$model == "ttest", ]$shift,
       appli_res[appli_res$model == "glmm", ]$shift,
       paired = TRUE) -> appli_ttest.vs.glmm

t.test(appli_res[appli_res$model == "ttest", ]$shift,
       appli_res[appli_res$model == "ehmos", ]$shift,
       paired = TRUE) -> appli_ttest.vs.ehmos
