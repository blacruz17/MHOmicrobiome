library(curatedMetagenomicData)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridtext)
library(janitor)
rm(list = ls())

dim(sampleMetadata)
length(unique(sampleMetadata$study_name))



## Westernized / Non-Westernized ###############################################

eur <- c("USA", "CAN",  
         "ITA", "GBR", "SWE", "DEU", "AUT", "DNK", "LUX", "FRA",
         "NOR", "SVK", "HUN", "EST", "FIN", "ISL", "IRL", "ESP",
         "NLD", "SLV", "RUS")
non.eur <- c("FJI", "CHN", "IDN", "MYS", "BRN", "KAZ", "BGD",
             "IND", "PHL", "MNG", "CMR", "PER", "MDG", "TZA",
             "LBR", "GHA", "ETH", "JPN", "KOR", "ISR", "SGP")
#check:
any(unique(sampleMetadata$country)[
      !unique(sampleMetadata$country) %in% c(eur, non.eur)
      ]) # should be false


dotplot.country <- sampleMetadata %>%
  select(study_name, country, non_westernized) %>%
  distinct() %>%
  mutate(eur = if_else(country %in% eur & non_westernized == "no", 1, 0)) %>%
  group_by(study_name) %>%
  summarise(eur = any(eur > 0), .groups = "drop") %>%
  arrange(eur)

# how many samples & studies does this exclude?:
sampleMetadata %>%
  mutate(eur = if_else(country %in% eur & non_westernized == "no", 1, 0)) %>%
  filter(eur == 0) %>%
  dim()

sum(dotplot.country$eur == FALSE)

## stool / no #########################
dotplot.stool <- sampleMetadata %>%
  select(study_name, body_site) %>%
  distinct() %>%
  group_by(study_name) %>%
  summarise(stool = any(body_site == "stool"), .groups = "drop") 


sampleMetadata %>%
  filter(body_site != "stool") %>%
  dim()


sum(dotplot.stool$stool == FALSE)

## antibiotics  #########################
dotplot.antb <- sampleMetadata %>%
  select(study_name, antibiotics_current_use) %>%
  distinct() %>%
  group_by(study_name) %>%
  summarise(antb = any(antibiotics_current_use == "no" | is.na(antibiotics_current_use)), 
            .groups = "drop",
            ) 

sampleMetadata %>%
  filter(antibiotics_current_use == "yes") %>%
  dim()

sum(dotplot.antb$antb == FALSE)


## age
dotplot.age <- sampleMetadata  %>%
  select(study_name, age) %>%
  group_by(study_name) %>%
  summarise(age = any(age >= 18 | is.na(age)), 
            .groups = "drop")
sum(dotplot.age$age == FALSE)

sampleMetadata %>%
  filter(age < 18) %>%
  dim()


## get all together and build dotplot
dotplot.data.1 <- dotplot.age %>%
  full_join(dotplot.antb, by = "study_name") %>%
  full_join(dotplot.country, by = "study_name") %>%
  full_join(dotplot.stool, by = "study_name") %>%
  mutate(total = age + antb + eur + stool) %>%
  arrange(desc(total))


table(dotplot.data.1$total)

order.dotplot <- dotplot.data.1 %>%
  mutate(order = case_when(
    total == 4 ~ 1,
    total == 3 & stool == 1 & antb == 1 & eur == 1 ~ 2, # solo 1
    total == 3 & stool == 1 & antb == 1 & eur == 0 ~ 3,
    total == 3 & stool == 1 & antb == 0 ~ 4,
    total == 3 & stool == 0 ~ 5,  # solo 1
    total == 2 & stool == 1 & antb == 1 ~ 6,
    total == 2 & stool == 1 & eur == 1 ~ 7,
    total == 2 & stool == 1 & age == 1 ~ 8,
    total == 2 & stool == 0 & eur == 0 ~ 9,
    total == 2 & stool == 0 & age == 0 ~ 10,
    total == 2 & stool == 0 & antb == 0 ~ 11,
    total == 1 & stool == 1 ~ 12,
    total == 1 & stool == 0 ~ 13
  )) %>%
  arrange(order) %>%
  pull(study_name)

panel.dotplot <- dotplot.data.1 %>%
  mutate(order = case_when(
    total == 4 ~ 1,
    total < 4 & stool == 1 ~ 2, 
    total < 4 & stool == 0 ~ 3)) %>%
  arrange(order) %>%
  pull(study_name)

dotplot.data.1 <- dotplot.data.1 %>%
  mutate(order = case_when(
    total == 4 ~ 1,
    total < 4 & stool == 1 ~ 2, 
    total < 4 & stool == 0 ~ 3)) %>%
  arrange(order)

dotplot.data <- dotplot.data.1 %>%
  select(study_name, stool, antb, age, eur, order) %>%
  reshape2::melt(id.vars = c("study_name", "order")) 

dotplot.data$study_name <- factor(dotplot.data$study_name,
                                     levels = rev(panel.dotplot))
dp <- dotplot.data %>% 
  ggplot(aes(y = study_name, x = variable, fill = value)) + 
  geom_point(shape = 21, color = "#264653") + 
  scale_fill_manual(values = c("white", "#264653"))

dp

dp.mx <- as.matrix(pivot_wider(dotplot.data,
                               id_cols = c("study_name", order), 
                               names_from = "variable", 
                               values_from = "value") %>% 
                     column_to_rownames(var = "study_name"))
dp.mx.2 <- apply(dp.mx, 2, as.numeric)
rownames(dp.mx.2) <- rownames(dp.mx)
dp.mx <- dp.mx.2
rm(dp.mx.2)

dp.mx <- dp.mx[order.dotplot, ]


## see all eligible samples ###############################
elig.df <- sampleMetadata %>% 
  # Westernized
  mutate(eur = if_else(country %in% eur & non_westernized == "no", 1, 0)) %>%
  filter(eur == 1) %>%
  # Stool samples
  filter(body_site == "stool" ) %>%
  # Age
  filter(age >= 18 | is.na(age)) %>%
  # Antibiotics
  filter(antibiotics_current_use == "no" | is.na(antibiotics_current_use)) %>%
  janitor::remove_empty()

dim(elig.df)

length(unique(elig.df$sample_id))
length(unique(elig.df$study_name))

# eligible samples w/ metadata: ###############################
elig.samples <- elig.df %>%
  select(sample_id, study_name,
         BMI, systolic_p, dyastolic_p,
         fasting_glucose, glucose, triglycerides, hdl,
         treatment) %>%
  reshape2::melt(id.vars = c("sample_id", "study_name")) %>%
  group_by(sample_id) %>%
  summarise(values = any(!is.na(value))) %>%
  filter(values == TRUE)

length(unique(elig.samples$sample_id))

# number of studies:
elig.df %>%
  filter(sample_id %in% elig.samples$sample_id) %>%
  pull(study_name) %>%
  unique() %>%
  length()


# keep only eligible samples:
elig.df <- elig.df %>%
  filter(sample_id %in% elig.samples$sample_id) 


# diseases: ##############
sort(unique(elig.df$study_condition))

non.elig.df <- elig.df %>% 
  filter(study_condition %in% 
           c("adenoma", "asthma", "CRC", "CDI", "IBD", "ME/CFS", 
             "melanoma", "migraine", "PD", "T1D") |
         study_name == "ChuDM_2017") %>% 
  remove_empty() 

dim(non.elig.df)

# number of studies remaining:
sort(unique(non.elig.df$study_name))


elig.df.2 <- elig.df %>% 
  filter(!study_condition %in% 
           c("adenoma", "asthma", "CRC", "CDI", "IBD", "ME/CFS", 
             "melanoma", "migraine", "PD", "T1D"),
         study_name != "ChuDM_2017") %>% 
  remove_empty() 

# how many studies were left out?
sum(!unique(non.elig.df$study_name) %in% unique(elig.df.2$study_name))
unique(non.elig.df$study_name)[!unique(non.elig.df$study_name) %in% unique(elig.df.2$study_name)]

unique(elig.df.2$study_name)
length(unique(elig.df.2$sample_id))

# "AsnicarF_2021" UNAVAILABLE 
# "HMP_2012" UNAVAILABLE
# "LeChatelierE_2013" UNAVAILABLE    
# "LifeLinesDeep_2016" UNAVAILABLE   
# "HMP_2019_ibdmdb" UNAVAILABLE          
# "NagySzakalD_2017"  UNAVAILABLE, chronic fatigue
unavail <- c("AsnicarF_2021", "HMP_2012",
             "LeChatelierE_2013", "NielsenHB_2014", 
             # Nielsen included here because it has samples in common with LeChatelier
             "LifeLinesDeep_2016",
             "HMP_2019_ibdmdb", "NagySzakalD_2017")
elig.df.2 %>%
  filter(study_name %in% unavail) %>%
  dim()

aa <- elig.df.2 %>%
  filter(study_name %in% unavail) %>%
  pull(sample_id) %>% unique() 
length(aa) # 2940 samples: measured but unavailable

bb <- elig.df.2 %>%
  filter(!study_name %in% unavail) %>%
  pull(sample_id) %>% unique() 

cc <- elig.df.2 %>%
  pull(sample_id) %>% unique() 

elig.df.3 <- elig.df.2 %>% filter(!study_name %in% unavail)
length(unique(elig.df.3$study_name))
length(unique(elig.df.3$sample_id)) # 6007 - 2940 = 3067 remaining

nonmeasured <- c("BedarfJR_2017", "DeFilippisF_2019", "FerrettiP_2018",
                 "HanniganGD_2017", "HansenLBS_2018", "IjazUZ_2017",
                 "KeohaneDM_2020", "LiJ_2014", 
                 "Obregon-TitoAJ_2015", "RaymondF_2016",
                 "SankaranarayananK_2015", "SchirmerM_2016", "ThomasAM_2018a",
                 "ThomasAM_2018b", "VogtmannE_2016", "WirbelJ_2018","XieH_2016",
                 "ZellerG_2014"  )
length(nonmeasured)

elig.df.3 %>%
  filter(study_name %in% nonmeasured) %>%
  dim()

elig.df.3 %>%
  filter(study_name %in% nonmeasured) %>%
  pull(sample_id) %>%
  unique() %>%
  length() # 1545 samples


# "BedarfJR_2017"  NON-MEASURED 
# "DeFilippisF_2019"  NON-MEASURED   # vegans!  
# "FerrettiP_2018"   NON-MEASURED
# "HanniganGD_2017"  NON-MEASURED     
# "HansenLBS_2018"  NON-MEASURED/ No healthy
# "IjazUZ_2017"  NON-MEASURED         
# "KeohaneDM_2020" NON-MEASURED
# "LiJ_2014"  No healthy, not enough metadata 
# "NielsenHB_2014"  Not found  
# "Obregon-TitoAJ_2015" NON-MEASURED  
# "RaymondF_2016" NON-MEASURED         
# "SankaranarayananK_2015" NON-MEASURED
# "SchirmerM_2016"  Not found
# "ThomasAM_2018a"  Not found       
# "ThomasAM_2018b"  Not found
# "VogtmannE_2016"  NON-MEASURED       
# "WirbelJ_2018"  NON-MEASURED         
# "XieH_2016"     NON-MEASURED
# "ZellerG_2014"  NON-MEASURED 

elig.df.4 <- elig.df.3 %>%
  filter(!study_name %in% nonmeasured)
dim(elig.df.4)
unique(elig.df.4$study_name)
length(unique(elig.df.4$sample_id)) 

# "Heitz-BuschartA_2016"  family, other diseases       
# "VincentC_2016" hospitalized patients        

elig.df.4 %>%
  filter(study_name %in% c("Heitz-BuschartA_2016",
                           "VincentC_2016")) %>%
  dim()

elig.df.4 %>%
  filter(study_name %in% c("Heitz-BuschartA_2016",
                           "VincentC_2016")) %>%
  pull(sample_id) %>% unique() %>% length()

final.df <- elig.df.4 %>%
  filter(!study_name %in% c("Heitz-BuschartA_2016",
                           "VincentC_2016"))

length(final.df$sample_id) 

final.df %>%
  select(study_name, sample_id) %>%
  distinct() %>%
  group_by(study_name) %>%
  count()


################################################################################
# HEATMAP ######################################################################
################################################################################

data.hm <- sampleMetadata %>%
  select(study_name,
         BMI, systolic_p, dyastolic_p, gender,
         fasting_glucose, glucose, triglycerides, hdl,
         treatment) %>%
 
  mutate(glucose_new = if_else(!is.na(fasting_glucose),
                               fasting_glucose, glucose),
         blood_pressure = if_else(!is.na(systolic_p),
                                  systolic_p, dyastolic_p)) %>%
  select(-c(fasting_glucose, glucose,
            systolic_p, dyastolic_p)) %>%
  reshape2::melt(id.vars = "study_name") %>%
  group_by(study_name, variable) %>%
  summarise(var_n = sum(!is.na(value)), .groups = "drop")


total_n <- sampleMetadata %>%
  select(study_name, sample_id) %>% count(study_name)

data.hm <- data.hm %>% 
  left_join(total_n, by = "study_name") %>%
  mutate(perc = (100 * var_n/n))


# Heatmap
hm.mx <- data.hm %>% 
  pivot_wider(id_cols = "study_name",
               names_from = "variable",
               values_from = "perc") %>%
  column_to_rownames(var = "study_name") %>%
  as.matrix()



################################################################################
# COMBINE BOTH PLOTS ###########################################################
################################################################################

colnames(dp.mx)[2:5] <- c("Stool samples", "Antibiotics", "Adults", "Location")

chosen_studies <- c("FengQ_2015", "HMP_2019_t2d", "KarlssonFH_2013", "MetaCardis_2020_a")
ha <- HeatmapAnnotation(Chosen = if_else(rownames(dp.mx) %in% chosen_studies, 
                                         "Yes", "No"),
                        col = list(Chosen = c("Yes" = "#7F3B08", 
                                             "No" = "#FEE0B6")),
                        which = "row",
                        gp = gpar(col = "gray60"),
                        name = "Chosen studies",
                        annotation_name_gp = gpar(fontface = "bold", 
                                                  fontsize = 10),
                        annotation_name_side = "top",
                        annotation_name_rot = 45,
                        annotation_legend_param =  list(title_gp = gpar(fontsize = 10,
                                                                        fontface = "bold"),
                                                        labels_gp = gpar(fontsize = 7)),
                        border = TRUE)


ht.1 <- Heatmap(dp.mx[ ,2:5],
                col = c(  "#FFF7F3", "#91017A"),
                rect_gp = gpar(col = "gray50", lwd = 1),
                border = TRUE,
                
                row_title = "Study", row_title_side = "left",
                row_title_gp = gpar(fontface = "bold",
                                    fontsize = 10),
                
                column_title = gt_render( "Inclusion <br> criteria"), 
                column_title_side = "bottom",
                column_title_gp = gpar(fontface = "bold",
                                       fontsize = 10),
                
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                row_split = factor(dp.mx[, 1], levels = 1:3),
                cluster_row_slices = FALSE,
                
                right_annotation = ha,
                
                cluster_columns = FALSE,
                
                column_names_side = "top",
                column_names_rot = 45,
                column_names_gp = gpar(fontsize = 10),
                
                show_row_names = FALSE,
                row_gap = unit(2, "mm"),
                
                heatmap_legend_param = list(title = "Elegibility",
                                            title_gp = gpar(fontsize = 10,
                                                            fontface = "bold"),
                                            labels_gp = gpar(fontsize = 7)),
                width = 1)

ht.1




pal <-  colorRampPalette(brewer.pal(n = 9, name = "RdPu"))(100)
colnames(hm.mx) <- c("BMI", "Gender", "Triglycerides", "HDL cholesterol",
                     "Medication", "Glucose", "Blood pressure")
hm.mx <- hm.mx[order.dotplot, c("BMI", "Triglycerides", "HDL cholesterol",
                                "Glucose", "Blood pressure", "Medication",
                                "Gender")]

ht.2 <-Heatmap(hm.mx,
        col = pal,
        rect_gp = gpar(col = "gray50", lwd = 1),
        border = TRUE,
        
        row_title = "Study", 
        row_title_side = "left",
        row_title_gp = gpar(fontface = "bold",
                            fontsize = 10),
        
        column_title = gt_render("Variables <br> of interest"),
        column_title_side = "bottom",
        column_title_gp = gpar(fontface = "bold",
                               fontsize = 10),
        
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        
        row_split = factor(dp.mx[, 1], levels = 1:3),
        cluster_row_slices = FALSE,
        
        column_names_side = "top",
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10),
        
        show_row_names = FALSE,
        
        heatmap_legend_param = list(title = "Availability",
                                    title_gp = gpar(fontsize = 10,
                                                    fontface = "bold"),
                                    labels_gp = gpar(fontsize = 7)),
        width = 6,
        row_gap = unit(2, "mm")
)
ht.2

ht_list <- ht.2 + ht.1
draw(ht_list, 
     ht_gap = unit(2, "mm"))

png("supp_heatmap.png",
    width = 15, height = 18, units = "cm", res=1200)
draw(ht_list, 
     ht_gap = unit(2, "mm"))
dev.off()