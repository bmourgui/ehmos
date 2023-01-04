##############
##
## Functions used to customize graph
##
## bastien.mourguiart@gmail.com
## file name: graph_functions.R
##
###############

library(ggplot2)

# Create a common theme for ggplot
ggtheme <- function(g){
  g +
    theme_bw() +
    theme(text = element_text(size=8) ,
          axis.title.x = element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = c("#1B9E77", "#D95F02"),
                      name = "Optimum type",
                      labels = c("Edge", "Middle")) +
    scale_color_manual(values = c("#1B9E77", "#D95F02"),
                       name = "Optimum type",
                       labels = c("Edge", "Middle")) +
    scale_shape_manual(values=c(25,21),
                       name = "Ecological type",
                       labels = c("Generalist", "Specialist")) +
    scale_x_discrete(labels=c("ehmos" = "EHMOS", "glmm" = "cGLMM",
                              "t.test" = "t-test")) -> g2
  return(g2)
}


# Nested plot with combinations if the 3 sub-scenarios

ggnested <- function(g, ylab = "RMSE (m)"){

  g +
    geom_violin(fill = "grey",
                position = position_dodge(0.9),
                alpha=0.5
                ) +
    ggh4x::facet_nested(marginality ~ sampling + specialization,
                        labeller = labeller(marginality = c("B1" = "B1: Mid sp only",
                                                            "B2" = "B2: Mid & Edge sp"),
                                            sampling = c("A1" = "A1: Uniform sampling",
                                                         "A2" = "A2: Unbalanced sampling"),
                                            specialization = c("C1" = "C1: Specialist only",
                                                               "C2" = "C2: Specialist & Generalist")
                                            )
                        ) +
    geom_jitter(aes(fill = T_opti,
                    color = T_opti,
                    shape = T_shape
                    ),
                size = 2,
                width = 0.2,
                alpha = 0.4
                ) +
      stat_summary(fun.data = "mean_se") +
      ylab(ylab) +
      xlab("")
}


