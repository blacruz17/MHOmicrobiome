library(tidyverse)
library(jsonlite)
library(gtsummary)
library(gt)
library(ggtext)

# K-Cores ######################################################################
muo <- read_csv("../data/MUO_kcores.csv") %>% mutate(class = "MUO")
mho <- read_csv("../data/MHO_kcores.csv") %>% mutate(class = "MHO")
muno <- read_csv("../data/MUNO_kcores.csv") %>% mutate(class = "MUNO")
mhno <- read_csv("../data/MHNO_kcores.csv") %>% mutate(class = "MHNO")

df <- rbind(muo, mho, muno, mhno) %>% 
  janitor::clean_names()


(df_total <- df %>%
  group_by(class) %>%
  summarise(total = n()))

(df_kcore <- df %>%
  group_by(class, k_core) %>%
  summarise(n = n()) %>%
  distinct() %>%
  arrange(class, k_core) %>%
  left_join(df_total, by = "class") %>%
    mutate(freq_kcore = 100 * n/total) %>%
    select(class, k_core, freq_kcore, total))

pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")

df_kcore  %>%
  ggplot(aes(x = k_core, y = freq_kcore)) +
  geom_col(aes(fill = class), col = "black") +
  scale_fill_manual(values = pal) +
  facet_wrap(. ~ class, 
             nrow = 1,
             strip.position = "top") +
  theme_minimal()  +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "none",
        axis.title = element_markdown()) +
  xlab("*k* values") +
  ylab("% nodes with maximum *k*") 


ggsave("../figures/kcores.png",
       width = 15,
       height = 5,
       unit = "cm",
       dpi = 1200)


# Node topology metrics ########################################################
shortest_paths <- fromJSON("../data/shortest_paths.json")
metrics <- read_csv("../data/network_metrics.csv")[, 2:6]

## Kruskal-Wallis --------------------------------------------------------------
metrics$red <- as.factor(metrics$Network)

kruskal_degree <- kruskal.test(degree ~ red, data = metrics)
kruskal_betweenness <- kruskal.test(betweenness ~ red, data = metrics)
kruskal_closeness <- kruskal.test(closeness ~ red, data = metrics)

print(kruskal_degree)
print(kruskal_betweenness)
print(kruskal_closeness)

# Shortest path lengths --------------------------------------------------------
df_paths <- do.call(rbind, lapply(names(shortest_paths), function(red) {
  data.frame(longitud = unlist(shortest_paths[[red]]), red = red)
}))

df_paths$red <- as.factor(df_paths$red)

kruskal_spl <- kruskal.test(longitud ~ red, data = df_paths)
print(kruskal_spl)

# p-adj ------------------------------------------------------------------------
p.adjust(c(kruskal_degree$p.value, kruskal_betweenness$p.value, 
           kruskal_closeness$p.value, kruskal_spl$p.value),
         method = "BH")