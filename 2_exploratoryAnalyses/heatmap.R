library(curatedMetagenomicData)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridtext)
rm(list = ls())


## Location

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

## stool / no
dotplot.stool <- sampleMetadata %>%
  select(study_name, body_site) %>%
  distinct() %>%
  group_by(study_name) %>%
  summarise(stool = any(body_site == "stool"), .groups = "drop") 


## antibiotics
dotplot.antb <- sampleMetadata %>%
  select(study_name, antibiotics_current_use) %>%
  distinct() %>%
  group_by(study_name) %>%
  summarise(antb = any(antibiotics_current_use == "no" | is.na(antibiotics_current_use)), 
            .groups = "drop",
            ) 


## age
dotplot.age <- sampleMetadata  %>%
  select(study_name, age) %>%
  group_by(study_name) %>%
  summarise(age = any(age >= 18 | is.na(age)), 
            .groups = "drop")


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

# Heatmap:

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


png("heatmap.png",
    width = 15, height = 18, units = "cm", res=1200)
draw(ht_list, 
     ht_gap = unit(2, "mm"))
dev.off()
