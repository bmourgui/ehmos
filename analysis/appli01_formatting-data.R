# R code for the application study in:
# A new method to explicitly estimate the shift of optimum along gradients in multispecies studies.
# B. Mourguiart, B. Liquet, K. Mengersen, T. Couturier, J. Mansons, Y. Braud, A. Besnard

# Script 1: Load species occurrence data and format them prior to model fitting


database <- readr::read_delim(here::here("data/BASE_PNM+Ventoux_2018_2019.csv"),
                            delim = "\t",
                            escape_double = FALSE,
                            trim_ws = TRUE)

data1980 <- readxl::read_excel(here::here("data", "BASE_data_PNM_gueguen.xlsx"),
                               sheet = "gueguen1980")

## 'Modern' data : From raw data to contingency table --------------------
database %>%
  dplyr::filter(stade != "juvenile") %>% #choice to omit juvenile as it is usually done and that s what probably did Gueguen
  dplyr::select(-c(operateur, rq_station, rq_ind, data_GPS, abondance, X_WGS84, Y_WGS84)) %>%
  #omit site with only two plots sampled
  dplyr::filter(station != "235(pas fini)") %>%
  #keep just the first fifth plots in sites 625 and A164 on which 10 plots have been sampled
  dplyr::filter(!(releve %in% c(540:544, 643:647))) %>%
  #omit phenological test conducted on the site A164 in 2018, keep only the sampling conducted the 9 aug 2018
  dplyr::filter(!(station == "A164" & date %in% c("23/07/2018", "17/09/2018"))) %>%
  dplyr::select(station, releve, altitude, taxon, vu, entendu, fauche)  -> data

data[is.na(data)] <- 0

data %>%
  dplyr::mutate(detection = as.numeric((entendu+vu+fauche)>0)) %>%
  dplyr::select(station, releve, taxon, detection) %>%
  dplyr::distinct() %>%  #one duplicated row, station 175 releve 621 taxon 'stenobothrus lineatus'
  tidyr::spread(key=taxon, value=detection) %>%
  dplyr::select(-c(`Gomphocerinae sp` , `Chorthippus sp`, ZERO_orthoptere)) -> tab #remove un-identified individuals, Gomphocerinae 2 and Chorthippus 1
tab[is.na(tab)] <- 0


## Homogenize taxon levels --------------------------------
# modern data
tab %>%
  dplyr::mutate("Anonconotus alpinus"=`Anonconotus baracunensis occidentalis`+ #Anonconotus alpinus was split in three taxa between 1980 and 2018
           `Anonconotus ghiliani`+`Anonconotus mercantouri`) %>%
  dplyr::select(-c(`Anonconotus baracunensis occidentalis`,
            `Anonconotus ghiliani`, `Anonconotus mercantouri`)) %>%
  dplyr::mutate("Calliptamus complex"=`Calliptamus italicus`+`Calliptamus siciliae`) %>% #high risk of missidentification between those 2 species, better treat them as a complex
  dplyr::select(-c(`Calliptamus italicus`, `Calliptamus siciliae`)) %>%
  dplyr::mutate("Chorthippus complex"= `Chorthippus biguttulus` + `Chorthippus brunneus brunneus` +
           `Chorthippus mollis` + `Chorthippus saulcyi daimei`) %>% #high risk of missidentification in this complex
  dplyr::select(-c(`Chorthippus biguttulus` , `Chorthippus brunneus brunneus` ,
            `Chorthippus mollis` , `Chorthippus saulcyi daimei`)) %>%
  dplyr::mutate("Podisma complex"= `Podisma dechambrei` + `Podisma pedestris`) %>% #high risk of missidentification between those 2 taxa in 1980, create a complex pedestris-dechambrei
  dplyr::select(-c(`Podisma dechambrei` , `Podisma pedestris`)) %>%
  dplyr::mutate("Yersinella complex"=`Yersinella beybienkoi`) %>%
  dplyr::select(-c(`Yersinella beybienkoi`)) -> data.modern

# historic data
data1980 %>%
  dplyr::mutate("Calliptamus complex"=`Calliptamus italicus`) %>% #high risk of missidentification with `Calliptamus siciliae`, better treat them as a complex
  dplyr::select(-c(`Calliptamus italicus`)) %>%
  dplyr::mutate("Chorthippus complex"= `Chorthippus biguttulus` +
           `Chorthippus mollis` + `Chorthippus saulcyi daimei`) %>% #high risk of missidentification in the complex biguttulus-brunneus-mollis-daimei", but brunneus not found or not identified as brunneus in 1980
  dplyr::select(-c(`Chorthippus biguttulus` ,
            `Chorthippus mollis` , `Chorthippus saulcyi daimei`)) %>%
  dplyr::mutate("Podisma complex"= `Podisma pedestris`) %>% #high risk of missidentification in the complex pedestris-dechambrei
  dplyr::select(-c(`Podisma pedestris`)) %>%
  dplyr::mutate("Yersinella complex"=`Yersinella raymondi`) %>% #possible confusion with Y. beybienkoi, gather the two species in one complex
  dplyr::select(-c(`Yersinella raymondi`)) -> data.historic

#View if there is incoherence in species names
matrix(c(colnames(data.historic),colnames(data.historic) %in% colnames(data.modern)),ncol=2) #also show that just two species had been found in 1980 and not in 2018-19, Barbitistes serricauda and Clonopsis gallica which was probably confounded with Pijnackeria masettii found in the same place in 2018
colnames(data.historic)[24] <- "Tessellana tessellata"
colnames(data.modern)[16] <- "Decticus verrucivorus"


## Merge datasets of occurrence data at the site-level ------------------------------
#transform datas in 'presence/absence' data at the site level, keep only shared species and re-sampled sites
data.modern %>%
  dplyr::group_by(station) %>%
  dplyr::summarise_if(is.numeric, .funs=function(x)(as.numeric(sum(x)>0))) %>%
  dplyr::select(intersect(colnames(data.modern), colnames(data.historic))) %>%
  dplyr::filter(station %in% data.historic$station) -> modern_data


data.historic %>%
  dplyr::select(-c(2:4)) %>%
  dplyr::mutate_at(.vars = -1, as.numeric) %>%
  dplyr::group_by(station) %>%
  dplyr::summarise_all(function(x)sum(x, na.rm=T)) %>%
  dplyr::mutate_if(is.numeric, .funs=function(x)(as.numeric(x>0))) %>%
  dplyr::select(intersect(colnames(data.modern), colnames(data.historic))) %>%
  dplyr::filter(station %in% data.modern$station) -> historic_data

data.frame("species"=colnames(historic_data[,-1]),
           "occurrence_hist"=colSums(historic_data[,-1]),
           "occurrence_mod"=colSums(modern_data[,-1]), row.names = 1:41) %>%
  dplyr::mutate("sufficient_occ"=(occurrence_hist>(0.05*135) & occurrence_mod>(0.05*135))) -> occu_tot

dplyr::bind_rows("historic" = historic_data, "modern" = modern_data, .id = "campaign") %>%
  #select species occurring in more than 5% of the sites in both period
  dplyr::select(c("campaign", "station", as.character(occu_tot[occu_tot$sufficient_occ == T,]$species))) -> dataF


## Reshape for modeling -----------------------
data %>%
  dplyr::select(station , releve, altitude) %>%
  dplyr::distinct(releve, .keep_all = T) %>%
  dplyr::group_by(station) %>%
  dplyr::summarise_at("altitude", mean) %>%
  dplyr::filter(station %in% dataF$station) %>%
  dplyr::inner_join(data.historic[,1:2], by="station") %>%
  dplyr::distinct() -> altitudes

dataF %>%
  tidyr::gather("species","Occ", -c(1:2)) %>%
  reshape2::melt(id.var=c("species", "station", "campaign"), measure.var="Occ") %>%
  reshape2::acast(station ~ species ~ campaign) -> Z

J <- dim(Z)[1] #number of sites
n <- dim(Z)[2] #number of species
S <- dim(Z)[3] #number of survey
alti <- altitudes$altitude.x
names(dimnames(Z)) <- list("site", "species", "period")
