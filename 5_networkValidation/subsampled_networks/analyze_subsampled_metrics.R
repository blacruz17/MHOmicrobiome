library(tidyverse)
library(RColorBrewer)
library(corrplot)
library(phyloseq)
library(gtsummary)
library(gt)
library(rstatix)
library(ggpubr)

rm(list = ls())

# 1. Data import ###############################################################
sampling.df <- read_csv("./randomSubsamples.csv")

physeq <- readRDS("./physeq_adj.rds")
metadata <- data.frame(sample_data(physeq)) %>%
  mutate(sex = if_else(sex == "Female", 1, 0)) %>%
  select(sex, age, bmi_kg_m2)

subsample.meta <- apply(sampling.df, 2, function(x) {
  women <- sum(metadata[x, "sex"])
  mean_bmi <- mean(metadata[x, "bmi_kg_m2"])
  mean_age <- mean(metadata[x, "age"])
  out <- c(women, mean_bmi, mean_age)
  names(out) <- c("Women", "BMI (Mean)", "Age (Mean)")
  return(out)
  }
  ) %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "netcode")

net.data <- read_csv("./subsampledNetsData.csv") %>%
  mutate(netcode = paste0(MetObesity, "_", number_network))

# 2. Get Stats #################################################################

tbl_todo <- net.data %>%
  select(-c(file, number_network, netcode)) %>%
  tbl_summary(by = "MetObesity",
              missing = "no") %>%
  add_n() %>%
  add_p() %>%
  add_q() %>%
  bold_p(q = TRUE)
tbl_todo


## ----pairwise--------------------------------------------------------------------------------------------------------
dunn.tb <-  net.data %>%
  select(-c(file, number_network, netcode))  %>% 
  reshape2::melt(id.vars = "MetObesity") %>% 
  group_by(variable) %>%
  dunn_test(value ~ MetObesity,
            p.adj = "BH")

dunn.tb

## ----wranglingPvals--------------------------------------------------------------------------------------------------
stat.test <- net.data %>%
  select(-c(file, number_network, netcode, max_KCore)) %>% 
  reshape2::melt(id.vars = "MetObesity") %>% 
  mutate(variable2 = variable) %>%
  group_by(variable2) %>%
  dunn_test(value ~ MetObesity,
            p.adj = "BH") %>%
  add_xy_position(scales = "free_y")


## plot ---------------
var.labs <- c("Order", "Size", "Edge density", "Clustering coefficient (mean)",
              "Shortest path length (mean)", "Betweenness centrality (mean)",
              "Closeness centrality (mean)", "Mean degree", "Maximum degree")
names(var.labs) <- unique(stat.test$variable2)


pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")

p <- net.data %>%
  select(-max_KCore) %>%
  reshape2::melt(id.vars = c("file", "MetObesity", 
                             "number_network", "netcode")) %>%
  rename(variable2 = variable) %>%
  ggplot(aes(x = MetObesity, y = value)) +
  geom_boxplot(aes(fill = MetObesity), 
               alpha=.6,
               outliers = FALSE,
               color = "black",
               # width = .1
  ) +
  geom_jitter(aes(shape = MetObesity, fill = MetObesity),
              width = .2,
              shape = 21,
              color = "gray10",
              alpha = .8,
              size = 2) +
  scale_color_manual(values = pal[c(1, 3, 4)],
                     aesthetics = c("fill", "color"),
                     name = "") +
  facet_wrap(. ~ variable2, 
             ncol = 3,
             scales = "free_y",
             labeller = labeller(variable2 = var.labs)
             ) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  stat_pvalue_manual(stat.test,
                     hide.ns = TRUE,
                     tip.length = 0.01,
                     bracket.shorten = 0.05,
                     step.increase = .1,
                     bracket.size = 0.3,
                     label.size = 2.5
  ) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    plot.title = element_text(face = 'bold'),
    axis.text = element_text(size = 12),
    strip.text = element_text(face = 'bold', size = 14,
                              color = "black"),
    strip.background = element_rect(
      color="white", fill="white",  linetype="solid"
    ))
p
ggsave("./subsampledNets_boxplots_stats.png",
       plot = p,
       bg = "white",
       height = 18,
       width = 20,
       units = "cm",
       dpi = 1200)
