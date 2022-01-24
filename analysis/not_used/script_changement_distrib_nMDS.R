#########################################;
##### Comparaison des distributions #####
#########################################;

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(Hmisc)
library(ggpubr)
library(stringr)
library(ggrepel)




load(file="data_rapport.RData")
gueguen1980 <- read_excel("data_PNM_gueguen.xlsx", sheet="gueguen1980")
gueguen2018 <- read_excel("data_PNM_gueguen.xlsx", sheet="2018")


gueguen1980 %>%
  filter(station %in% gueguen2018$station)  %>% #filtre les stations communes aux deux periodes
  select(-c(1,4:5)) %>% #garde seulement le nom des stations, l'altitude et les esp
  arrange(station) -> gueguen1980


gueguen1980[,-c(1,2)] -> tab.mds1980 #enleve nom et altitude
tab.mds1980[tab.mds1980>1] <- 1 #passe de la densite a pres/abs

gueguen2018  %>%
  select(-1) -> tab.mds2018 


esp.communes <- intersect(colnames(tab.mds1980), colnames(tab.mds2018)) 


tab.mds1980 %>%
  select(esp.communes) -> tab.mds1980
tab.mds2018 %>%
  select(esp.communes) -> tab.mds2018

tab.mds2018 <- tab.mds2018[rowSums(tab.mds2018)>0,] #possible que en ne conservant que les esp communes, des lignes soient nulles. On les enleve pour pouvoir faire l'analyse multivar
tab.mds1980 <- tab.mds1980[rowSums(tab.mds1980)>0,]


#### Approche multivariee ; nMDS 2D + lm ####

#nMDS sur la matrice de distance de jaccard
set.seed(188)
mds1980 <- metaMDS(tab.mds1980, k=2, trymax=250, dist = "jaccard") 
mds2018 <- metaMDS(tab.mds2018, k=2, trymax=250, dist = "jaccard", shrink=FALSE) 


NMDS.data2018 <- as.data.frame(mds2018$species[colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4,]) #seules les especes non rares sont representees
NMDS.data1980 <- as.data.frame(mds1980$species[colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4,])

as.data.frame(stringr::str_split_fixed(colnames(tab.mds2018), " ",3)) %>%
  select(-3) %>%
  mutate("taxon"=paste(str_sub(V1,1,3),str_sub(V2,1,3), sep="_")) %>%
  select(taxon) -> sp.names

NMDS.data1980[,3] <- sp.names[colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4,]
NMDS.data2018[,3] <- sp.names[colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4,]
colnames(NMDS.data1980)[3] <- "taxon"
colnames(NMDS.data2018)[3] <- "taxon"

#Relation (lineaire) positionnement graphique et altitude
points.mds1980 <- data.frame("altitude"=gueguen1980$altitude , #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
                             "MDS1"=mds1980$points[,1], "MDS2"=mds1980$points[,2])
points.mds2018 <- data.frame("altitude"=cov.sites$altitude , #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
                             "MDS1"=mds2018$points[,1], "MDS2"=mds2018$points[,2])
lm1980 <- lm(altitude ~ MDS1 + MDS2, data=points.mds1980)
lm2018 <- lm(altitude ~ MDS1 + MDS2, data=points.mds2018)

#Prediction position altitudinale des especes
sp.predict1980 <- predict.lm(lm1980, newdata = NMDS.data1980, interval = "confidence")
sp.predict2018 <- predict.lm(lm2018, newdata = NMDS.data2018, interval = "confidence")

data.frame("periode"=as.factor(rep(c("1983-1988",2018), each=nrow(NMDS.data1980))),
           "taxon"=rep(NMDS.data2018$taxon,2), 
           "altitude"= c(sp.predict1980[,"fit"],sp.predict2018[,"fit"]),
           "altitude.min"= c(sp.predict1980[,"lwr"],sp.predict2018[,"lwr"]),
           "altitude.max"= c(sp.predict1980[,"upr"],sp.predict2018[,"upr"]),
           "methode"="multivar") -> lm.altiSP

lm.altiSP %>% 
  ggplot(aes(x=taxon, y=altitude, color=`periode`)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=altitude.min, ymax=altitude.max), width=0.5) +
  theme_classic() +
  ylab("Altitude (m)") +
  theme(axis.text.x = element_text(face = "italic", angle = -45, vjust = 0.5, hjust=0, size=14),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14)) +
  scale_color_manual(values = c("cyan","darkred"))



#### Approche multivariee ; nMDS 3D + lm ####

#nMDS sur la matrice de distance de jaccard
set.seed(188)
mds1980_3D <- metaMDS(tab.mds1980, k=3, trymax=250, dist = "jaccard") 
mds2018_3D <- metaMDS(tab.mds2018, k=3, trymax=250, dist = "jaccard", shrink=FALSE) 


NMDS.data2018_3D <- as.data.frame(mds2018_3D$species[colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4,]) #seules les especes non rares sont representees
NMDS.data1980_3D <- as.data.frame(mds1980_3D$species[colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4,])

NMDS.data1980_3D[,4] <- sp.names[colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4,]
NMDS.data2018_3D[,4] <- sp.names[colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4,]
colnames(NMDS.data1980_3D)[4] <- "taxon"
colnames(NMDS.data2018_3D)[4] <- "taxon"

#Relation (lineaire) positionnement graphique et altitude
points.mds1980_3D <- data.frame("altitude"=gueguen1980$altitude , #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
                             "MDS1"=mds1980_3D$points[,1], "MDS2"=mds1980_3D$points[,2], "MDS3"=mds1980_3D$points[,3])
points.mds2018_3D <- data.frame("altitude"=cov.sites$altitude , #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
                             "MDS1"=mds2018_3D$points[,1], "MDS2"=mds2018_3D$points[,2], "MDS3"=mds1980_3D$points[,3])
lm1980_3D <- lm(altitude ~ MDS1 + MDS2 + MDS3, data=points.mds1980_3D)
lm2018_3D <- lm(altitude ~ MDS1 + MDS2 + MDS3, data=points.mds2018_3D)

#Prediction position altitudinale des especes
sp.predict1980_3D <- predict.lm(lm1980_3D, newdata = NMDS.data1980_3D, interval = "confidence")
sp.predict2018_3D <- predict.lm(lm2018_3D, newdata = NMDS.data2018_3D, interval = "confidence")



data.frame("periode"=as.factor(rep(c("1983-1988",2018), each=nrow(NMDS.data1980_3D))),
           "taxon"=rep(NMDS.data2018_3D$taxon,2), 
           "altitude"= c(sp.predict1980_3D[,"fit"],sp.predict2018_3D[,"fit"]),
           "altitude.min"= c(sp.predict1980_3D[,"lwr"],sp.predict2018_3D[,"lwr"]),
           "altitude.max"= c(sp.predict1980_3D[,"upr"],sp.predict2018_3D[,"upr"]),
           "methode"="multivar_3D") -> lm.altiSP_3D

lm.altiSP_3D %>% 
  ggplot(aes(x=taxon, y=altitude, color=`periode`)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=altitude.min, ymax=altitude.max), width=0.5) +
  theme_classic() +
  ylab("Altitude (m)") +
  theme(axis.text.x = element_text(face = "italic", angle = -45, vjust = 0.5, hjust=0, size=14),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14)) +
  scale_color_manual(values = c("cyan","darkred"))



#### Approche comparaison de moyennes ####
tab.alti.sp2018 <- tab.mds2018[,colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4]*cov.sites$altitude #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
tab.alti.sp1980 <- tab.mds1980[,colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4]*cov.sites$altitude #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici

tab.alti.sp1980[,colSums(tab.alti.sp1980)>4 & colSums(tab.alti.sp2018)>4]


tab.alti.sp1980[tab.alti.sp1980==0] <- NA
tab.alti.sp2018[tab.alti.sp2018==0] <- NA




mean.alti.sp2018 <- apply(tab.alti.sp2018,2,function(x)mean(x, na.rm=TRUE))
mean.alti.sp1980 <- apply(tab.alti.sp1980,2,function(x)mean(x, na.rm=TRUE))

confint.mean <- function(x){
  n <- sum(!is.na(x))
  borne.inf <- mean(x, na.rm=TRUE) - 1.96*(sd(x, na.rm=TRUE)/sqrt(n))
  borne.sup <- mean(x, na.rm=TRUE) + 1.96*(sd(x, na.rm=TRUE)/sqrt(n))
  int <- c(borne.inf, borne.sup)
  return(int)
}

confint.alti.sp2018 <- apply(tab.alti.sp2018,2,confint.mean)
confint.alti.sp1980 <- apply(tab.alti.sp1980,2,confint.mean)

mean.altiSP <- data.frame("periode"=as.factor(rep(c("1983-1988",2018), each=nrow(sp.names))),
                          "taxon"=sp.names,
                          "altitude"=c(mean.alti.sp1980, mean.alti.sp2018), 
                          "altitude.min"=c(confint.alti.sp1980[1,], confint.alti.sp2018[1,]),
                          "altitude.max"=c(confint.alti.sp1980[2,], confint.alti.sp2018[2,]),
                          "methode"="moyenne")

mean.altiSP %>%
  filter(taxon %in% lm.altiSP$taxon) -> mean.altiSP

ggplot(mean.altiSP, aes(x=taxon, y=altitude, color=`periode`)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=altitude.min, ymax=altitude.max), width=0.5) +
  theme_classic() +
  ylab("Altitude (m)") +
  theme(axis.text.x = element_text(face = "italic", angle = -45, vjust = 0.5, hjust=0, size=14),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14)) +
  scale_color_manual(values = c("cyan","darkred"))


rbind(mean.altiSP, lm.altiSP, lm.altiSP_3D) %>%
  ggplot(aes(x=taxon, y=altitude, color=`periode`)) +
  facet_wrap(.~methode) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=altitude.min, ymax=altitude.max), width=0.5) +
  theme_classic() +
  ylab("Altitude (m)") +
  theme(axis.text.x = element_text(face = "italic", angle = -45, vjust = 0.5, hjust=0, size=14),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14)) +
  scale_color_manual(values = c("cyan","darkred")) +
  theme(panel.grid.major = element_line(linetype = c("28")))





##### Approche un seul tableau #####

tab.mdsUNIQUE <- cbind(tab.mds1980, tab.mds2018)
colnames(tab.mdsUNIQUE) <- c(paste(sp.names$taxon, "1980", sep=""),paste(sp.names$taxon, "2018", sep=""))

mds.unique <- metaMDS(tab.mdsUNIQUE, k=2, trymax=250, dist = "jaccard")
NMDS.dataunique <- as.data.frame(mds.unique$species[colSums(tab.mds1980)>4 & colSums(tab.mds2018)>4,]) #seules les especes non rares sont representees


NMDS.dataunique[,3] <- rownames(NMDS.dataunique)

colnames(NMDS.dataunique)[3] <- "taxon"

#Relation (lineaire) positionnement graphique et altitude
points.mdsUNIQUE <- data.frame("altitude"=gueguen1980$altitude , #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
                                "MDS1"=mds.unique$points[,1], "MDS2"=mds.unique$points[,2])


plot.mdsUNIQUE <- ggplot(NMDS.dataunique, aes(y = MDS1, x = MDS2)) +
  geom_point(size = 4, alpha = 0.8) + 
  geom_point(data=as.data.frame(mds.unique$points), aes(y = MDS1, x = MDS2, color=cov.sites$altitude), 
             size=4, shape=17) +
  geom_text_repel(aes(label=taxon), size=6) + #plots the NMDS points, with shape by topo type
  theme_bw() + #for aesthetics
  labs(colour = "Altitude") + #another way to set the labels, in this case, for the colour legend
  scale_colour_gradient(high = "#A91101", low = "yellow") + #here we set the high and low of the colour scale.  Can delete to go back to the standard blue, or specify others
  theme(legend.key = element_blank(),  #removes the box around each legend item
        legend.position = "bottom", #legend at the bottom
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.just = "centre") +
  scale_fill_manual(values = c("white","red","blue","green","black")) +
  ylim(c(-2.1,2))

plot.mdsUNIQUE


data.frame("taxon"=NMDS.dataunique[1:27,3], 
      "MDS1.1980"=NMDS.dataunique[1:27,1],
      "MDS1.2018"=NMDS.dataunique[28:54,1], 
      "delta.MDS1"=NMDS.dataunique[28:54,1]-NMDS.dataunique[1:27,1],
      "delta.alti"=lm.altiSP$altitude[28:54]-lm.altiSP$altitude[1:27],
      "delta.moy"=mean.altiSP$altitude[28:54]-mean.altiSP$altitude[1:27]) %>%
  arrange(-delta.MDS1)








#################################
########## SIMULATIONS ##########
#################################


###### Basic simulation #####

# 6 species 
#     spA specialist of low altitudes, upslope shift of 300m
#     spB specialist of high altitudes, downslope shift of 300m
#     spC specialist of medium altitudes, no shift
#     spD generalist, upslope shift of 300m
#     spE generalist, downslope shift of 300m
#     spF generalist, no shift

###### 77 sites on one transect each 25m of elevation ######

alti <- seq(900, 2800, 25)
station <- 1:77

spA.t1 <- c(rbinom(9,1,0.90), rep(0,68))
spA.t2 <- c(rep(0,12), rbinom(9,1,0.90), rep(0,56))

spB.t1 <- c(rep(0,68), rbinom(9,1,0.90))
spB.t2 <- c(rep(0,56), rbinom(9,1,0.90), rep(0,12))

spC.t1 <- c(rep(0,34), rbinom(9,1,0.90), rep(0,34))
spC.t2 <- c(rep(0,34), rbinom(9,1,0.90), rep(0,34))

spD.t1 <- c(rbinom(47,1,0.90), rep(0,30))
spD.t2 <- c(rep(0,12), rbinom(47,1,0.90), rep(0,18))

spE.t1 <- c(rep(0,30), rbinom(47,1,0.90))
spE.t2 <- c(rep(0,18), rbinom(47,1,0.90), rep(0,12))

spF.t1 <- c(rep(0,15), rbinom(47,1,0.90), rep(0,15))
spF.t2 <- c(rep(0,15), rbinom(47,1,0.90), rep(0,15))

X.t1 <- cbind(spA.t1, spB.t1, spC.t1, spD.t1, spE.t1, spF.t1)
X.t2 <- cbind(spA.t2, spB.t2, spC.t2, spD.t2, spE.t2, spF.t2)


##### mean comparison #####

X.t1b <- X.t1*alti
X.t2b <- X.t2*alti

X.t1b[X.t1b==0] <- NA
X.t2b[X.t2b==0] <- NA

mean.alti.t1 <- apply(X.t1b,2,function(x)mean(x, na.rm=TRUE))
mean.alti.t2 <- apply(X.t2b,2,function(x)mean(x, na.rm=TRUE))

confint.mean <- function(x){
  n <- sum(!is.na(x))
  borne.inf <- mean(x, na.rm=TRUE) - 1.96*(sd(x, na.rm=TRUE)/sqrt(n))
  borne.sup <- mean(x, na.rm=TRUE) + 1.96*(sd(x, na.rm=TRUE)/sqrt(n))
  int <- c(borne.inf, borne.sup)
  return(int)
}

confint.alti.t1 <- apply(X.t1b,2,confint.mean)
confint.alti.t2 <- apply(X.t2b,2,confint.mean)

mean.altiSP <- data.frame("periode"=as.factor(rep(c("t1","t2"), each=ncol(X.t1))),
                          "taxon"=colnames(X.t1),
                          "altitude"=c(mean.alti.t1, mean.alti.t2), 
                          "altitude.min"=c(confint.alti.t1[1,], confint.alti.t2[1,]),
                          "altitude.max"=c(confint.alti.t1[2,], confint.alti.t2[2,]),
                          "methode"="moyenne")


ggplot(mean.altiSP, aes(x=taxon, y=altitude, color=`periode`)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=altitude.min, ymax=altitude.max), width=0.5) +
  theme_classic() +
  ylab("Altitude (m)") +
  theme(axis.text.x = element_text(face = "italic", angle = -45, vjust = 0.5, hjust=0, size=14),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14)) +
  scale_color_manual(values = c("cyan","darkred"))



##### (mds2D + lm) method #####

set.seed(188)

X.t1MDS <- X.t1[rowSums(X.t1)>0,]
X.t2MDS <- X.t2[rowSums(X.t2)>0,]

mds.t1 <- metaMDS(X.t1MDS, k=2, trymax=250, dist = "jaccard") 
mds.t2 <- metaMDS(X.t2MDS, k=2, trymax=250, dist = "jaccard") 


NMDS.datat1 <- as.data.frame(mds.t1$species) 
NMDS.datat2 <- as.data.frame(mds.t2$species)

NMDS.datat1[,3] <- colnames(X.t1)
NMDS.datat2[,3] <- colnames(X.t1)

colnames(NMDS.datat1)[3] <- "taxon"
colnames(NMDS.datat2)[3] <- "taxon"

#Relation (lineaire) positionnement graphique et altitude
points.mdst1 <- data.frame("altitude"=alti[rowSums(X.t1)>0] , #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
                             "MDS1"=mds.t1$points[,1], "MDS2"=mds.t1$points[,2])
points.mdst2 <- data.frame("altitude"=alti[rowSums(X.t2)>0] , #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
                             "MDS1"=mds.t2$points[,1], "MDS2"=mds.t2$points[,2])
lmt1 <- lm(altitude ~ MDS1 + MDS2, data=points.mdst1)
lmt2 <- lm(altitude ~ MDS1 + MDS2, data=points.mdst2)

#Prediction position altitudinale des especes
sp.predictt1 <- predict.lm(lmt1, newdata = NMDS.datat1, interval = "confidence")
sp.predictt2 <- predict.lm(lmt2, newdata = NMDS.datat2, interval = "confidence")

data.frame("periode"=as.factor(rep(c("t1","t2"), each=nrow(NMDS.datat1))),
           "taxon"=rep(NMDS.datat2$taxon,2), 
           "altitude"= c(sp.predictt1[,"fit"],sp.predictt2[,"fit"]),
           "altitude.min"= c(sp.predictt1[,"lwr"],sp.predictt2[,"lwr"]),
           "altitude.max"= c(sp.predictt1[,"upr"],sp.predictt2[,"upr"]),
           "methode"="multivar") -> lm.altiSP

lm.altiSP %>% 
  ggplot(aes(x=taxon, y=altitude, color=`periode`)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=altitude.min, ymax=altitude.max), width=0.5) +
  theme_classic() +
  ylab("Altitude (m)") +
  theme(axis.text.x = element_text(face = "italic", angle = -45, vjust = 0.5, hjust=0, size=14),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14)) +
  scale_color_manual(values = c("cyan","darkred"))


##### same sites as in the PNM #####


sp1.tA <- (cov.sites$altitude > 1200) & (cov.sites$altitude < 1800)
sp2.tA <- (cov.sites$altitude > 1800) & (cov.sites$altitude < 2500)
sp3.tA <- (cov.sites$altitude > 1100) & (cov.sites$altitude < 1900)
sp4.tA <- (cov.sites$altitude > 1900) & (cov.sites$altitude < 2600)

sp1.tB <- (cov.sites$altitude > 1500) & (cov.sites$altitude < 2100)
sp2.tB <- (cov.sites$altitude > 1500) & (cov.sites$altitude < 2200)
sp3.tB <- (cov.sites$altitude > 1800) & (cov.sites$altitude < 2100)
sp4.tB <- (cov.sites$altitude > 2000) & (cov.sites$altitude < 2700)

X.tA <- apply(cbind(sp1.tA, sp2.tA, sp3.tA, sp4.tA),2,as.numeric)
X.tB <- apply(cbind(sp1.tB, sp2.tB, sp3.tB, sp4.tB),2,as.numeric)


##### mean comparison #####

X.tAb <- X.tA*cov.sites$altitude
X.tBb <- X.tB*cov.sites$altitude

X.tAb[X.tAb==0] <- NA
X.tBb[X.tBb==0] <- NA

mean.alti.tA <- apply(X.tAb,2,function(x)mean(x, na.rm=TRUE))
mean.alti.tB <- apply(X.tBb,2,function(x)mean(x, na.rm=TRUE))

confint.mean <- function(x){
  n <- sum(!is.na(x))
  borne.inf <- mean(x, na.rm=TRUE) - 1.96*(sd(x, na.rm=TRUE)/sqrt(n))
  borne.sup <- mean(x, na.rm=TRUE) + 1.96*(sd(x, na.rm=TRUE)/sqrt(n))
  int <- c(borne.inf, borne.sup)
  return(int)
}

confint.alti.tA <- apply(X.tAb,2,confint.mean)
confint.alti.tB <- apply(X.tBb,2,confint.mean)

mean.altiSP <- data.frame("periode"=as.factor(rep(c("tA","tB"), each=ncol(X.tA))),
                          "taxon"=colnames(X.tA),
                          "altitude"=c(mean.alti.tA, mean.alti.tB), 
                          "altitude.min"=c(confint.alti.tA[1,], confint.alti.tB[1,]),
                          "altitude.max"=c(confint.alti.tA[2,], confint.alti.tB[2,]),
                          "methode"="moyenne")


ggplot(mean.altiSP, aes(x=taxon, y=altitude, color=`periode`)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=altitude.min, ymax=altitude.max), width=0.5) +
  theme_classic() +
  ylab("Altitude (m)") +
  theme(axis.text.x = element_text(face = "italic", angle = -45, vjust = 0.5, hjust=0, size=14),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14)) +
  scale_color_manual(values = c("cyan","darkred"))



##### (mds2D + lm) method #####


X.tAMDS <- X.tA[rowSums(X.tA)>0,]
X.tBMDS <- X.tB[rowSums(X.tB)>0,]

mds.tA <- metaMDS(X.tAMDS, k=2, trymax=250, dist = "jaccard") 
mds.tB <- metaMDS(X.tBMDS, k=2, trymax=250, dist = "jaccard") 


NMDS.datatA <- as.data.frame(mds.tA$species) 
NMDS.datatB <- as.data.frame(mds.tB$species)

NMDS.datatA[,3] <- colnames(X.tA)
NMDS.datatB[,3] <- colnames(X.tA)

colnames(NMDS.datatA)[3] <- "taxon"
colnames(NMDS.datatB)[3] <- "taxon"

#Relation (lineaire) positionnement graphique et altitude
points.mdstA <- data.frame("altitude"=alti[rowSums(X.tA)>0] , #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
                           "MDS1"=mds.tA$points[,1], "MDS2"=mds.tA$points[,2])
points.mdstB <- data.frame("altitude"=alti[rowSums(X.tB)>0] , #si on a enleve des lignes n ayant que de 0 pour les esp, penser a enelver les stations ici
                           "MDS1"=mds.tB$points[,1], "MDS2"=mds.tB$points[,2])
lmtA <- lm(altitude ~ MDS1 + MDS2, data=points.mdstA)
lmtB <- lm(altitude ~ MDS1 + MDS2, data=points.mdstB)

#Prediction position altitudinale des especes
sp.predicttA <- predict.lm(lmtA, newdata = NMDS.datatA, interval = "confidence")
sp.predicttB <- predict.lm(lmtB, newdata = NMDS.datatB, interval = "confidence")

data.frame("periode"=as.factor(rep(c("tA","tB"), each=nrow(NMDS.datatA))),
           "taxon"=rep(NMDS.datatB$taxon,2), 
           "altitude"= c(sp.predicttA[,"fit"],sp.predicttB[,"fit"]),
           "altitude.min"= c(sp.predicttA[,"lwr"],sp.predicttB[,"lwr"]),
           "altitude.max"= c(sp.predicttA[,"upr"],sp.predicttB[,"upr"]),
           "methode"="multivar") -> lm.altiSP

lm.altiSP %>% 
  ggplot(aes(x=taxon, y=altitude, color=`periode`)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=altitude.min, ymax=altitude.max), width=0.5) +
  theme_classic() +
  ylab("Altitude (m)") +
  theme(axis.text.x = element_text(face = "italic", angle = -45, vjust = 0.5, hjust=0, size=14),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14)) +
  scale_color_manual(values = c("cyan","darkred"))


