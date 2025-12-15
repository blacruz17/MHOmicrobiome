library(dplyr)
library(tidyr)
library(UpSetR)
library(vegan)
library(ggplot2)
library(magick)
library(RColorBrewer)
library(readr)
library(tibble)
library(factoextra)
library(FactoMineR)

rm(list = ls())


# Filenames #########################################################
fn_upset  <- "panel_upset_ncbi.png"
fn_pca   <- "panel_nmds_pca.png"
fn_final_png <- "Figure_Keystone_Composite.png"
df_upset <- "upsetdf_ncbi.csv"

key.data <- read_csv("keystone_data_subsampledNets.csv")

# Create matrix with all info ##################################################
key.mx <- key.data %>%
  select(Species, keystone, Network) %>%
  mutate(keystone = as.numeric(keystone)) %>%
  pivot_wider(
    names_from = "Network",
    values_from = "keystone",
    values_fill = 0
  ) %>%
  column_to_rownames(var = "Species") %>%
  as.matrix()

key.mx <- key.mx[rowSums(key.mx) != 0, ]

# --- Image settings ---
group.colors <- c(
  "MHNO" = "#264653",
  "MUO"  = "#edafb8",
  "MUNO" = "#703d67"
)

width_in  <- 6
height_in <- 5
res_dpi   <- 600

######################## 1) UpSet plot #########################################

upset.df <- key.data %>%
  select(Species, Group, keystone) %>%
  group_by(Species, Group) %>%
  summarise(keystone = as.numeric(max(keystone)), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = keystone, values_fill = 0)

upset.df <- as.data.frame(upset.df)
upset.df[, -1] <- lapply(upset.df[, -1], as.integer)

write_csv(upset.df, df_upset)

png(fn_upset,
    width  = width_in * res_dpi,
    height = height_in * res_dpi,
    res    = res_dpi)
par(mar = c(4, 4, 4, 2))

UpSetR::upset(
  upset.df,
  sets = c("MHNO", "MUNO", "MUO"),
  keep.order = TRUE,
  main.bar.color = "#1B4F72",
  sets.bar.color = "grey40",
  order.by = "freq",
  text.scale = c(1.8, 1.6, 1.4, 1.4, 1.6, 1.8)
)

dev.off()

######################## 2) PCA plot ###########################################
comm.mx <- t(key.mx)
res.pca <- PCA(comm.mx, graph = FALSE)

ind_df <- as.data.frame(res.pca$ind$coord)
ind_df$Sample <- rownames(ind_df)
ind_df$Group <- gsub("_[0-9]+", "", ind_df$Sample)

var1 <- round(res.pca$eig[1, 2], 1)
var2 <- round(res.pca$eig[2, 2], 1)

p_pca <- ggplot(ind_df, aes(x = Dim.1, y = Dim.2)) +
  stat_ellipse(aes(color = Group, fill = Group),
               type = "t", level = 0.95, 
               linetype = 2,
               geom = "polygon",
               alpha = 0.2
  ) +
  geom_point(aes(fill = Group),
             shape = 21, size = 3, color = "black") +
  scale_fill_manual(values = group.colors) +
  scale_color_manual(values = group.colors) +
  labs(
    x = paste0("PC1 (", var1, "%)"),
    y = paste0("PC2 (", var2, "%)"),
    fill = "Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    text = element_text(family = "sans")
  )
p_pca

ggsave(fn_pca, p_pca, width = width_in, height = height_in, dpi = res_dpi)

######################## 3) Composite figure ###################################

img_up   <- image_read(fn_upset)
img_nmds <- image_read(fn_pca)

target_height_px <- height_in * res_dpi
img_nmds <- image_scale(img_nmds, paste0("x", target_height_px))
img_up   <- image_scale(img_up,  paste0("x", target_height_px))

img_up   <- image_annotate(img_up, "B", size = 150, weight = 700, 
                           color = "black", location = "+40+80")
img_nmds <- image_annotate(img_nmds, "A", size = 150, weight = 700, 
                           color = "black", location = "+40+80")

combined <- image_append(c(img_nmds, img_up), stack = FALSE)

image_write(combined, path = fn_final_png, format = "png",
            density = res_dpi, quality = 100)