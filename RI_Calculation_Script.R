####### Script to calculate RI #######

# Installation of R package to calculate linear retention index (RI).

# Installation of "MetaboCoreUtils" package
#install.packages("remotes")
#remotes::install_github("rformassspectrometry/MetaboCoreUtils")

# Loading "MetaboCoreUtils" library
library("MetaboCoreUtils")

# Extra libraries
library("readxl")
library("writexl")

# Experimental RI for features that match with NIST libraries
## Loadding the retention time (RT) of each n-alkane
alkane_rt <- c(3.4665, 5.0805, 7.2105, 9.6470, 12.1855, 14.7160, 17.1685,
               19.5140, 21.7500, 23.8885, 25.9295, 28.1660, 31.0990, 35.0805,
               39.1355, 41.8615, 43.9890, 45.7755, 47.3495, 48.7680, 50.0745,
               51.2945, 52.4410, 53.6010, 54.9075)
## Loadding the RI of each n-alkane
alkane_ri <- c(900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
               1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800,
               2900, 3000, 3100, 3200, 3300)
## n-alkanes DataFrame
alkane_data <- data.frame(rtime = alkane_rt, rindex = alkane_ri)


## Experimental Retention Index (RI) calculation
################################# MS-DIAL ######################################
##Feature list deconvoluted with MS-DIAL
msdial_rt <- read_excel("../B_grandiflora_metabolomics/Data/B_gradiflora_MSDIAL_Feat_list.xlsx", sheet = 1)
##RI calculation
msdial_ri <- indexRtime(msdial_rt$`Retention Time`, alkane_data)
##Adding the RI to the MS-DIAL feature list
msdial_rt$RI_using_R <- msdial_ri
## Moving the RI column close to the RT column
msdial_rt <- msdial_rt[,c(1:6, 13, 7:12)]
##Exporting the MS-DIAL feature list with RI
write_xlsx(msdial_rt, "../B_grandiflora_metabolomics/Result/B_gradiflora_MSDIAL_Feat_list_RI.xlsx")


################################# MZmine #######################################
##Feature list deconvoluted with MZmine (3 to 43 min dataset)
mzmine_rt <- read_excel("../B_grandiflora_metabolomics/Data/B_gradiflora_MZmine_Feat_list_3to43min.xlsx", sheet = 1)
##RI calculation (3 to 43 min dataset)
mzmine_ri <- indexRtime(mzmine_rt$`row retention time`, alkane_data)
##Adding the RI to the MZmine feature list (3 to 43 min dataset)
mzmine_rt$RI <- mzmine_ri 
## Moving the RI column close to the RT column (3 to 43 min dataset)
mzmine_rt <- mzmine_rt[,c(1:2, 47, 3:46)]
##Exporting the MZmine feature list with RI
write_xlsx(mzmine_rt, "../B_grandiflora_metabolomics/Result/B_gradiflora_MZmine_Feat_list_3to43min_RI.xlsx")

## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with MZmine (43 to 53 min dataset)
mzmine_rt43to53 <- read_excel("../B_grandiflora_metabolomics/Data/B_gradiflora_MZmine_Feat_list_43to53min.xlsx", sheet = 1)
##RI calculation (43 to 53 min dataset)
mzmine_ri43to53 <- indexRtime(mzmine_rt43to53$`row retention time`, alkane_data)
##Adding the RI to the MZmine feature list (43 to 53 min dataset)
mzmine_rt43to53$RI <- mzmine_ri43to53 
## Moving the RI column close to the RT column (43 to 53 min dataset)
mzmine_rt43to53 <- mzmine_rt43to53[,c(1:3, 48, 4:47)]
##Exporting the MZmine feature list with RI
write_xlsx(mzmine_rt43to53, "../B_grandiflora_metabolomics/Result/B_gradiflora_MZmine_Feat_list_43to53min_RI.xlsx")


################################# eRah #########################################
##Feature deconvoluted with eRah
rtime_erah <- c(time...
                )
indexRtime(rtime_erah, alkane_data)


################################# MSHub ########################################
##Feature deconvoluted with MSHub
rtime_mshub <- c(time...
                 )
indexRtime(rtime_mshub, alkane_data)

