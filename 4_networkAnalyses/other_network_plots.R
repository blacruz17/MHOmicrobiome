library(SpiecEasi)
library(phyloseq)
library(igraph)
library(patchwork)
library(tidyverse)
library(viridis)
library(rstatix)
library(ggpubr)
library(ggtext)
library(gtsummary)

rm(list = ls())

# Calculate NR50/pc differences (Random networks) ---------------------------------
random.nr50 <- read_csv("./random_nr50.csv")

random.nr50.long <- reshape2::melt(random.nr50)
kruskal.test(random.nr50.long$value, random.nr50.long$variable) # significant!
colnames(random.nr50.long) <- c("MetObesity", "NR50")

pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")


stat.test <- random.nr50.long %>%
  dunn_test(NR50 ~  MetObesity,
            p.adj = "BH") %>%
  add_y_position()

pNR50 <- random.nr50.long %>%
  ggplot(aes(x = MetObesity, y = NR50)) +
  geom_boxplot(aes(fill = MetObesity),
               outliers = FALSE,
               alpha = .5) +
  geom_jitter(aes(color = MetObesity),
              width = .1) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  stat_pvalue_manual(stat.test,
                     hide.ns = TRUE,
                     tip.length = 0.01,
                     bracket.shorten = 0.05,
                     step.increase = .01,
                     bracket.size = 0.3,
                     label.size = 4
  ) +
  theme_bw()


# Percolation threshold:
random.fc <- read_csv("./random_fc.csv")
random.fc.long <- reshape2::melt(random.fc)
kruskal.test(random.fc.long$value, random.fc.long$variable) # significant!
colnames(random.fc.long) <- c("MetObesity", "fc")

pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")


stat.test <- random.fc.long %>%
  dunn_test(fc ~  MetObesity,
            p.adj = "BH") %>%
  add_y_position()

pFc <- random.fc.long %>%
  ggplot(aes(x = MetObesity, y = fc)) +
  geom_boxplot(aes(fill = MetObesity),
               outliers = FALSE,
               alpha = .5) +
  geom_jitter(aes(color = MetObesity),
              width = .1) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  stat_pvalue_manual(stat.test,
                     hide.ns = TRUE,
                     tip.length = 0.01,
                     bracket.shorten = 0.05,
                     step.increase = .01,
                     bracket.size = 0.3,
                     label.size = 4
  ) +
  theme_bw()


pBoth <- ggarrange(pNR50, pFc, common.legend = TRUE, legend = FALSE)
pBoth

ggsave("./figRandom_NR50.png",
       plot = pBoth,
       width = 15,
       height = 5,
       unit = "cm",
       dpi = 1200)

ggsave("./figRandom_fc.png",
       plot = random.fc.long %>%
         ggplot(aes(x = MetObesity, y = fc)) +
         geom_boxplot(aes(fill = MetObesity),
                      outliers = FALSE,
                      alpha = .5) +
         geom_jitter(aes(color = MetObesity),
                     width = .1, 
                     size = .1) +
         scale_fill_manual(values = pal) +
         scale_color_manual(values = pal) +
         theme_bw() + 
         labs(x = "Network", y = "Percolation threshold") + 
         theme(legend.position = "none"),
       width = 7,
       height = 5,
       unit = "cm",
       dpi = 600)



# K-cores ----------------------------------------------------------------------

kcores.mhno <- read_csv("./MHNO_kcores.csv")
kcores.mho <- read_csv("./MHO_kcores.csv")
kcores.muno <- read_csv("./MUNO_kcores.csv")
kcores.muo <- read_csv("./MUO_kcores.csv")

kcore.df <- bind_rows(
  kcores.mhno %>% mutate(Network = "MHNO"),
  kcores.mho %>% mutate(Network = "MHO"),
  kcores.muno %>% mutate(Network = "MUNO"),
  kcores.muo %>% mutate(Network = "MUO"))

colnames(kcore.df) <- c("Node", "KCore", "Network")

kcore.df %>%
  mutate(KCore = as.factor(KCore)) %>%
  count(Network, KCore) %>%
  group_by(Network) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(percentage = 100 * n/total)


kcore.df %>%
  mutate(KCore = as.factor(KCore)) %>%
  count(Network, KCore) %>%
  group_by(Network) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(percentage = 100 * n/n)  %>%
  ggplot(aes(x = KCore, y = n)) +
  geom_col(aes(fill = Network), col = "black") +
  scale_fill_manual(values = pal) +
  facet_wrap(. ~ Network, 
             nrow = 1,
             strip.position = "top") +
  labs(
    x = "*k* values",
    y = "% nodes with <br> maximum *k*"
  ) +
  theme_classic()  +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "none",
        axis.title = element_markdown(size = 10),
        axis.text = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )


ggsave("./kcores_20251128.png",
       width = 10,
       height = 5,
       unit = "cm",
       dpi = 1200)


## stats ---------------------------------------------------------------------
shortest_paths <- fromJSON("../results/shortest_paths.json")
metrics <- read_csv("../results/network_metrics.csv")[, 2:6]

metrics$red <- as.factor(metrics$Network)

kruskal_degree <- kruskal.test(degree ~ red, data = metrics)
kruskal_betweenness <- kruskal.test(betweenness ~ red, data = metrics)
kruskal_closeness <- kruskal.test(closeness ~ red, data = metrics)

print(kruskal_degree)
print(kruskal_betweenness)
print(kruskal_closeness)

df_paths <- do.call(rbind, lapply(names(shortest_paths), function(red) {
  data.frame(longitud = unlist(shortest_paths[[red]]), red = red)
}))

df_paths$red <- as.factor(df_paths$red)

kruskal_spl <- kruskal.test(longitud ~ red, data = df_paths)
print(kruskal_spl)

p.adjust(c(kruskal_degree$p.value, kruskal_betweenness$p.value, 
           kruskal_closeness$p.value, kruskal_spl$p.value),
         method = "BH")

# Network subsampling analyses -------------------------------------------------
net100 <- read_csv("./subsampledNetsData.csv") 

## Stats -----------------------------------------------------------------------
tbl_todo <- net100  %>%
  select(-max_KCore, -max_degree, -file, -number_network) %>%
  tbl_summary(by = "MetObesity",
              missing = "no") %>%
  add_n() %>%
  add_p() %>%
  add_q() %>%
  bold_p(q = TRUE)
tbl_todo
# all significantly different

# pairwise Dunn tests:
dunn.tb <- net100  %>%
  select(-max_KCore, -max_degree, -file, -number_network) %>% 
  reshape2::melt(id.vars = "MetObesity") %>% 
  group_by(variable) %>%
  dunn_test(value ~ MetObesity,
            p.adj = "BH")

dunn.tb

# get p-value table for stats_pvalue_manual:
stat.test <- net100  %>%
  select(-max_KCore, -max_degree, -file, -number_network) %>% 
  reshape2::melt(id.vars = "MetObesity") %>% 
  mutate(variable2 = variable) %>%
  group_by(variable2) %>%
  dunn_test(value ~ MetObesity,
            p.adj = "BH") %>%
  add_xy_position(scales = "free_y")

var.labs <- c("Number CCs", "Order", "Size", "% Positive Edges",
              "Edge Density", "ACC", "ASPL", "Mean Betweenness",
              "Mean Closeness", "Mean Degree")

names(var.labs) <- unique(stat.test$variable2)

## Plot ------------------------------------------------------------------------
p.net100 <- net100 %>%
  select(-max_KCore, -max_degree) %>%
  reshape2::melt(id.vars = c("file", "MetObesity", "number_network")) %>%
  rename(variable2 = variable) %>%
  ggplot(aes(x = MetObesity, y = value)) +
  geom_boxplot(aes(fill = MetObesity),
               outliers = FALSE,
               alpha = .5) +
  geom_jitter(aes(color = MetObesity),
              width = .1,
              size = .5) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  facet_wrap(. ~ variable2, 
             ncol = 3,
             scales = "free_y",
             labeller = labeller(variable2 = var.labs)) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  stat_pvalue_manual(stat.test,
                     hide.ns = TRUE,
                     step.increase = 0.05,
                     tip.length = 0.01,
                     bracket.shorten = 0.05,
                     bracket.size = 0.3,
                     label.size = 2.5
  ) +
  labs(x = "Network", y = "Metric values") +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    plot.title = element_text(face = 'bold'),
    axis.text = element_text(size = 6),
    strip.text = element_text(face = 'bold', size = 7,
                              color = "black"),
    strip.background = element_rect(
      color="white", fill="white",  linetype="solid"
    ))
p.net100

ggsave("./fig100plots.png",
       plot = p.net100,
       width = 15,
       height = 15,
       unit = "cm",
       dpi = 1200)

