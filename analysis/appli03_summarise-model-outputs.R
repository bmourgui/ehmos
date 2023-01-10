# R code for the application study in:
# A new method to explicitly estimate the shift of optimum along gradients in multispecies studies.
# B. Mourguiart, B. Liquet, K. Mengersen, T. Couturier, J. Mansons, Y. Braud, A. Besnard

# Script 3: Transform models' outputs into summary statistics (e.g., mean shift estimates)

# (Outputs of models used in the published manuscript can be download at: https://doi.org/10.5061/dryad.0gb5mkm59)


## Extract mcmc samples of ecological parameters
ehmos <- mcmc_ehmos(out = out_appli_ehmos,
                    scale = TRUE,
                    x = alti)
glmm <- mcmc_glmm(out = out_appli_glmm,
                    scale = TRUE,
                    x = alti)
## Derive point estimates and CrI of species-specific ecological parameters
summary_glmm <- summary_mcmc(X = glmm)
summary_ehmos <- summary_mcmc(X = ehmos)
summary_ttest <- list(
  "pmax" = NULL,
  "width" = t.test_width(Z, alti),
  "opt" = apply(Z[,,1], 2, function(x)(mean.method(x, alti))),
  "shift" = apply(Z, 2, function(x)(mean.method2(x, alti)))
)

## Create summary data frame
appli_res <- data.frame("model" = rep(c("ehmos", "glmm", "ttest"), each = n),
                        "sp" = dimnames(Z)$species,
                        "n1" = apply(Z[ , , 1], 2, sum),
                        "n2" = apply(Z[ , , 2], 2, sum),

                        "opt" = c(summary_ehmos$opt[2, ], summary_glmm$opt[2, ], summary_ttest$opt[2, ]),
                        "lwr.opt" = c(summary_ehmos$opt[1, ], summary_glmm$opt[1, ], summary_ttest$opt[1, ]),
                        "upr.opt" = c(summary_ehmos$opt[3, ], summary_glmm$opt[3, ], summary_ttest$opt[3, ]),

                        "shift" = c(summary_ehmos$shift[2, ], summary_glmm$shift[2, ], summary_ttest$shift[2, ]),
                        "lwr.shift" = c(summary_ehmos$shift[1, ], summary_glmm$shift[1, ], summary_ttest$shift[1, ]),
                        "upr.shift" = c(summary_ehmos$shift[3, ], summary_glmm$shift[3, ], summary_ttest$shift[3, ]),

                        "width1" = c(summary_ehmos$width[2, , 1], summary_glmm$width[2, , 1], summary_ttest$width[, 1]), # just an observation with ttest method
                        "lwr.width1" = c(summary_ehmos$width[1, , 1], summary_glmm$width[1, , 1], rep(NA, n)),
                        "upr.width1" = c(summary_ehmos$width[3, , 1], summary_glmm$width[3, , 1],  rep(NA, n)),

                        "width2" = c(summary_ehmos$width[2, , 2], summary_glmm$width[2, , 2], summary_ttest$width[ , 2]),
                        "lwr.width2" = c(summary_ehmos$width[1, , 2], summary_glmm$width[1, , 2],  rep(NA, n)),
                        "upr.width2" = c(summary_ehmos$width[3, , 2], summary_glmm$width[3, , 2],  rep(NA, n)),

                        "pmax1" = c(summary_ehmos$pmax[2, , 1], summary_glmm$pmax[2, , 1], rep(NA, n)),
                        "lwr.pmax1" = c(summary_ehmos$pmax[1, , 1], summary_glmm$pmax[1, , 1], rep(NA, n)),
                        "upr.pmax1" = c(summary_ehmos$pmax[3, , 1], summary_glmm$pmax[3, , 1], rep(NA, n)),

                        "pmax2" = c(summary_ehmos$pmax[2, , 2], summary_glmm$pmax[2, , 2], rep(NA, n)),
                        "lwr.pmax2" = c(summary_ehmos$pmax[1, , 2], summary_glmm$pmax[1, , 2], rep(NA, n)),
                        "upr.pmax2" = c(summary_ehmos$pmax[3, , 2], summary_glmm$pmax[3, , 2], rep(NA, n))

                        )
