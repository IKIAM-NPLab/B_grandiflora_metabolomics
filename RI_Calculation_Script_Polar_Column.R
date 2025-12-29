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
alkane_rt <- c(3.855, 5.395, 7.323, 9.435, 11.623, 13.823, 15.973, 18.058,
               20.068, 22.001, 23.858, 25.645, 27.521, 29.863, 32.900, 36.808,
               39.885, 42.153, 44.221, 46.620, 49.536)
## Loadding the RI of each n-alkane
alkane_ri <- c(1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
               2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200)
## n-alkanes DataFrame
alkane_data <- data.frame(rtime = alkane_rt, rindex = alkane_ri)


## Experimental Retention Index (RI) calculation
################################# MS-DIAL ######################################
##Feature list deconvoluted with MS-DIAL (3 to 21 min dataset)
msdial_rt <- read_excel("../B_grandiflora_metabolomics/Data/MSDIAL_Feature_List_3to21min_Polar_Column.xlsx", sheet = 1)
##RI calculation (3 to 21 min dataset)
msdial_ri <- indexRtime(msdial_rt$`Average Rt(min)`, alkane_data)
##Adding the RI to the MS-DIAL feature list (3 to 21 min dataset)
msdial_rt$RI_using_R <- msdial_ri
## Moving the RI column close to the RT column (3 to 21 min dataset)
msdial_rt <- msdial_rt[,c(1:3, 69, 4:68)]
##Exporting the MS-DIAL feature list with RI (3 to 21 min dataset)
write_xlsx(msdial_rt, "../B_grandiflora_metabolomics/Result/MSDIAL_Feature_List_3to21min_Polar_Column_RI.xlsx")

##Feature list deconvoluted with MS-DIAL (21 to 37 min dataset)
msdial_rt21to37 <- read_excel("../B_grandiflora_metabolomics/Data/MSDIAL_Feature_List_21to37min_Polar_Column.xlsx", sheet = 1)
##RI calculation (21 to 37 min dataset)
msdial_ri21to37 <- indexRtime(msdial_rt21to37$`Average Rt(min)`, alkane_data)
##Adding the RI to the MS-DIAL feature list (21 to 37 min dataset)
msdial_rt21to37$RI_using_R <- msdial_ri21to37
## Moving the RI column close to the RT column (21 to 37 min dataset)
msdial_rt21to37 <- msdial_rt21to37[,c(1:3, 69, 4:68)]
##Exporting the MS-DIAL feature list with RI (21 to 37 min dataset)
write_xlsx(msdial_rt21to37, "../B_grandiflora_metabolomics/Result/MSDIAL_Feature_List_21to37min_Polar_Column_RI.xlsx")

##Feature list deconvoluted with MS-DIAL (37 to 44 min dataset)
msdial_rt37to44 <- read_excel("../B_grandiflora_metabolomics/Data/MSDIAL_Feature_List_37to44min_Polar_Column.xlsx", sheet = 1)
##RI calculation (37 to 44 min dataset)
msdial_ri37to44 <- indexRtime(msdial_rt37to44$`Average Rt(min)`, alkane_data)
##Adding the RI to the MS-DIAL feature list (37 to 44 min dataset)
msdial_rt37to44$RI_using_R <- msdial_ri37to44
## Moving the RI column close to the RT column (37 to 44 min dataset)
msdial_rt37to44 <- msdial_rt37to44[,c(1:3, 69, 4:68)]
##Exporting the MS-DIAL feature list with RI (37 to 44 min dataset)
write_xlsx(msdial_rt37to44, "../B_grandiflora_metabolomics/Result/MSDIAL_Feature_List_37to44min_Polar_Column_RI.xlsx")

################################# MZmine #######################################
##Feature list deconvoluted with MZmine (3 to 21 min dataset)
mzmine_rt <- read_excel("../B_grandiflora_metabolomics/Data/MZmine_Feature_List_3to21min_polar_column.xlsx", sheet = 1)
##RI calculation (3 to 21 min dataset)
mzmine_ri <- indexRtime(mzmine_rt$`row retention time`, alkane_data)
##Adding the RI to the MZmine feature list (3 to 21 min dataset)
mzmine_rt$RI <- mzmine_ri
## Moving the RI column close to the RT column (3 to 21 min dataset)
mzmine_rt <- mzmine_rt[,c(1:3, 47, 4:46)]
##Exporting the MZmine feature list with RI
write_xlsx(mzmine_rt, "../B_grandiflora_metabolomics/Result/MZmine_Feature_List_3to21min_polar_column_RI.xlsx")

##Feature list deconvoluted with MZmine (21 to 37 min dataset)
mzmine_rt21to37 <- read_excel("../B_grandiflora_metabolomics/Data/MZmine_Feature_List_21to37min_polar_column.xlsx", sheet = 1)
##RI calculation (21 to 37 min dataset)
mzmine_ri21to37 <- indexRtime(mzmine_rt21to37$`row retention time`, alkane_data)
##Adding the RI to the MZmine feature list (21 to 37 min dataset)
mzmine_rt21to37$RI <- mzmine_ri21to37
## Moving the RI column close to the RT column (21 to 37 min dataset)
mzmine_rt21to37 <- mzmine_rt21to37[,c(1:3, 48, 4:47)]
##Exporting the MZmine feature list with RI
write_xlsx(mzmine_rt21to37, "../B_grandiflora_metabolomics/Result/MZmine_Feature_List_21to37min_polar_column_RI.xlsx")

##Feature list deconvoluted with MZmine (37 to 44 min dataset)
mzmine_rt37to44 <- read_excel("../B_grandiflora_metabolomics/Data/MZmine_Feature_List_37to44min_polar_column.xlsx", sheet = 1)
##RI calculation (37 to 44 min dataset)
mzmine_ri37to44 <- indexRtime(mzmine_rt37to44$`row retention time`, alkane_data)
##Adding the RI to the MZmine feature list (37 to 44 min dataset)
mzmine_rt37to44$RI <- mzmine_ri37to44
## Moving the RI column close to the RT column (37 to 44 min dataset)
mzmine_rt37to44 <- mzmine_rt37to44[,c(1:3, 48, 4:47)]
##Exporting the MZmine feature list with RI
write_xlsx(mzmine_rt37to44, "../B_grandiflora_metabolomics/Result/MZmine_Feature_List_37to44min_polar_column_RI.xlsx")






################################# eRah #########################################
## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with eRah (3 to 21.10 min dataset)
erah_rt3to21 <- read_excel("../B_grandiflora_metabolomics/Data/B_gradiflora_eRah_Feat_list_3to21min_Polar_Column.xlsx", sheet = 1)
##RI calculation (3 to 21 min dataset)
erah_ri3to21 <- indexRtime(erah_rt3to21$tmean, alkane_data)
##Adding the RI to the eRah feature list (3 to 21 min dataset)
erah_rt3to21$RI <- erah_ri3to21
## Moving the RI column close to the RT column (3 to 21 min dataset)
erah_rt3to21 <- erah_rt3to21[,c(1:4, 40, 5:39)]
##Exporting the eRah feature list with RI
write_xlsx(erah_rt3to21, "../B_grandiflora_metabolomics/Result/B_gradiflora_eRah_Feat_list_3to21min_Polar_Column_RI.xlsx")






## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with eRah (21 to 37 min dataset)
erah_rt21to37 <- read_excel("../B_grandiflora_metabolomics/Data/B_gradiflora_eRah_Feat_list_21to37min_Polar_Column.xlsx", sheet = 1)
##RI calculation (21 to 37 min dataset)
erah_ri21to37 <- indexRtime(erah_rt21to37$tmean, alkane_data)
##Adding the RI to the eRah feature list (21 to 37 min dataset)
erah_rt21to37$RI <- erah_ri21to37
## Moving the RI column close to the RT column (21 to 37 min dataset)
erah_rt21to37 <- erah_rt21to37[,c(1:4, 40, 5:39)]
##Exporting the eRah feature list with RI
write_xlsx(erah_rt21to37, "../B_grandiflora_metabolomics/Result/B_gradiflora_eRah_Feat_list_21to37min_Polar_Column_RI.xlsx")







## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with eRah (37 to 44 min dataset)
erah_rt37to44 <- read_excel("E:/B_grandiflora/B_grandiflora_metabolomics/Data/B_gradiflora_eRah_Feat_list_37to44_Polar_Column.xlsx", sheet = 1)
##RI calculation (37 to 44 min dataset)
erah_ri37to44 <- indexRtime(erah_rt37to44$tmean, alkane_data)
##Adding the RI to the eRah feature list (37to 44 min dataset)
erah_rt37to44$RI <- erah_ri37to44
## Moving the RI column close to the RT column (37 to 44 min dataset)
erah_rt37to44 <- erah_rt37to44[,c(1:4, 40, 5:39)]
##Exporting the eRah feature list with RI
write_xlsx(erah_rt37to44, "E:/B_grandiflora/B_grandiflora_metabolomics/Result/B_gradiflora_eRah_Feat_list_37to44_Polar_Column_RI.xlsx")


















################################# MSHub ########################################
## Experimental Retention Index (RI) calculation
##Feature list deconvoluted with MShub
mshub_rt <- read_excel("E:/B_grandiflora/B_grandiflora_metabolomics/Data/B_gradiflora_MSHub_Feat_list.xlsx", sheet = 1)
##RI calculation
mshub_ri <- indexRtime(mshub_rt$`row retention time`, alkane_data)
##Adding the RI to the MSHub feature list
mshub_rt$RI <- mshub_ri
## Moving the RI column close to the RT column
mshub_rt <- mshub_rt[,c(1:4, 39, 5:38)]
##Exporting the eRah feature list with RI
write_xlsx(mshub_rt, "E:/B_grandiflora/B_grandiflora_metabolomics/Result/B_gradiflora_MSHub_Feat_list_RI.xlsx")

