###############################%
#
# Extended result analysis of the application
# Not used to produce 1st submitted MS
#
# file name: appli04b_extended-results.R
# file dependencies: appli03*.R and previous
#
# 29/12/2021
#
###############################%

source(here::here("R/analysis_functions.R"))
load(here::here("results/appli_out-ehmos.RData"))
load(here::here("results/appli_out-glmm.RData"))

#### 1. Graph of SRCs along sampling gradient ####
## 1.1. Create SRC for each species and model ----
occu_appli_ehmos <- pred.ehmos(out = out_appli_ehmos,
                               scaled = TRUE, # Was the environmental vector scaled when model ran? yes
                               x.pred = alti, # Env. vector on which occupancy probabilities will be predicted. Here, sampled altitudes
                               full_matrix = FALSE) # Do we wan full occupancy matrix? No, just summary (medians and CrI)

occu_appli_glmm <- pred.glmm(out = out_appli_glmm,
                               scaled = TRUE, # Was the environmental vector scaled when model ran? yes
                               x.pred = alti, # Env. vector on which occupancy probabilities will be predicted. Here, sampled altitudes
                               full_matrix = FALSE) # Do we wan full occupancy matrix? No, just summary (medians and CrI)


## 1.2. Create objects used in graphs ----
# Pres/Abs data to plot
Z.graph <- cbind(reshape2::melt(Z),
                 "alti" = alti) %>%
  dplyr::mutate("species" = rep(rep(1:24, each = 135),2))

# Species names to plot
Z.graph %>%
  tidyr::separate(species,
                  c("genus", "sp")) %>%
  dplyr::mutate("genus" = substr(genus, 1, 2),
                "sp.names" = paste(genus, sp, sep = ". ")) %>%
  dplyr::select(sp.names) %>%
  dplyr::distinct() -> sp.names

# Optima to plot
appli_res %>%
  dplyr::mutate("species" = rep(1:24,3),
                "2" = opt + shift,
                "1" = opt) %>%
  dplyr::select(species, model, `1`, `2`) %>%
  tidyr::gather("sampling_occ", "opt", 3:4) -> d.opt

## 1.3. Create graphs ----
ggplot() +
  geom_line(data = occu_appli_ehmos,
            aes(x = alti, y = `50%`),
            color = "#1B9E77") +
  geom_line(data = occu_appli_glmm,
            aes(x = alti, y = `50%`),
            color = "#7570B3") +
  geom_jitter(data = Z.graph,
              aes(x = alti, y = value),
              color = "grey",
              width = 0,
              height = 0.01,
              alpha = 0.5) +
  geom_point(data = d.opt,
             aes(x = opt,
                 y = 1.01,
                 color = model,
                 fill = model
                 ),
             size=3) +
  ylab("Occupancy probability") +
  xlab("Altitude") +
  facet_grid(sampling_occ ~ species)




#### 2. Inspect residuals ####

plot(out_appli_glmm$sims.list$p.fit ~ out_appli_glmm$sims.list$p.fitnew, xlim=c(800,940), ylim=c(800,940))
abline(0,1)

plot(out_appli_ehmos$sims.list$p.fit ~ out_appli_ehmos$sims.list$p.fitnew, xlim=c(800,940), ylim=c(800,940))
abline(0,1)

p.value_glmm <- 0
for(s in 1:length(out_appli_glmm$sims.list$p.fit)){
  if(out_appli_glmm$sims.list$p.fit[s]<out_appli_glmm$sims.list$p.fitnew[s]){
    p.value_glmm <- p.value_glmm + 1
  }
}
p.value_glmm <- p.value_glmm/length(out_appli_glmm$sims.list$p.fit)

p.value_ehmos <- 0
for(s in 1:length(out_appli_ehmos$sims.list$p.fit)){
  if(out_appli_ehmos$sims.list$p.fit[s]<out_appli_ehmos$sims.list$p.fitnew[s]){
    p.value_ehmos <- p.value_ehmos + 1
  }
}
p.value_ehmos <- p.value_ehmos/length(out_appli_ehmos$sims.list$p.fit)


p.value_glmmSP <- NULL
for (i in 1:24){
  p.value_glmmSP[i] <- 0
  for(s in 1:length(out_appli_glmm$sims.list$p.fit)){
    if(out_appli_glmm$sims.list$p.fitSP[s,i]<out_appli_glmm$sims.list$p.fitnewSP[s,i]){
      p.value_glmmSP[i] <- p.value_glmmSP[i] + 1
    }
  }
  print(p.value_glmmSP[i]/nrow(out_appli_glmm$sims.list$p.fitSP))
}


p.value_ehmosSP <- NULL
for (i in 1:24){
  p.value_ehmosSP[i] <- 0
  for(s in 1:length(out_appli_ehmos$sims.list$p.fit)){
    if(out_appli_ehmos$sims.list$p.fitSP[s,i]<out_appli_ehmos$sims.list$p.fitnewSP[s,i]){
      p.value_ehmosSP[i] <- p.value_ehmosSP[i] + 1
    }
  }
  print(p.value_ehmosSP[i]/nrow(out_appli_ehmos$sims.list$p.fitSP))
}
