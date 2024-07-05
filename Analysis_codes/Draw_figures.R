#################################################################################################
########## Assembly and maturation of calf gut microbiome from neonate to post-puberty ##########
#################################################################################################
# Code written by Jae-Yun Lee
# Revised at 2024-07-02

library(ggpubr)
library(ggtext)
library(ape)
library(cowplot)
library(phyloseq)
library(tidyverse)

# Load files --------------------------------------------------------------

table_path = "ASV_table.tsv"
tax_path = "Taxonomy.tsv"
tree_path = "Phylogenetic_tree.nwk"
md_path = "Metadata.tsv"

table <- read_delim(table_path) %>% column_to_rownames("ASV")
table <- otu_table(table, taxa_are_rows = TRUE)

tax.df <- read_delim(tax_path)
tax.df <- tax.df %>% column_to_rownames("Feature.ID")
tax.df <- tax_table(as.matrix(tax.df))

phy.tree <- read.tree(tree_path)

md <- read.table(md_path, sep = "\t", header = TRUE)
rownames(md) <- md$sample_name
md <- sample_data(md)
ps <- phyloseq(table, md, tax.df, phy.tree)

# Fig.2c-e: Depth plot of samples ------------------------------------------------------------

# Generating depth dataframe
depth.df <- data.frame(depth = sample_sums(ps)) # Calcaulte depth of each sample
depth.df$sample_name <- rownames(depth.df)
depth.df <- data.frame(md) %>% left_join(depth.df, by = "sample_name")

rpl.text <- c("Large_intestine" = "Li", "Rumen" = "Ru", "Small_intestine" = "Si")
depth.df <- depth.df %>% mutate(sample_origin = str_replace_all(sample_origin, rpl.text))

type.lev <- c("Feces", "Li", "Si", "Ru", "Milk", "Diet")
depth.df$sample_origin <- factor(depth.df$sample_origin, levels = type.lev)

sex.lev <- c("M", "F")
day.level <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9")
sample.level <- c("Feces", "Li", "Si", "Ru",  "Milk", "Diet")
depth.df$host_sex <- factor(depth.df$host_sex, levels = sex.lev)
depth.df$sample_origin <- factor(depth.df$sample_origin, levels = sample.level)

sex.col <- c(M = "blue", `F` = "red")
tp.cols <- c("#000000", "#004000", "#008000", "#00C000", "#01FE00", "#41BE00", 
             "#817E00", "#C13E00", "#FF6E00")
names(tp.cols) <- day.level
sample.cols <- c(Feces = "#307F00", Li = "#FFBF00", 
                 Si = "darkorange", Ru = "#FF0000", Milk = "#0F52BA", Diet = "violet")
depth.df$timepoints <- factor(depth.df$timepoints, levels = day.level)

# Fig. 2c
ggplot(depth.df, aes(x = sample_origin, y = log10(depth), color = sample_origin)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 21, fill = NA, width = 0.1) +
  theme_minimal_hgrid() +
  labs(x = NULL, color = "Sample", y = "log~10~Depth") +
  guides(color = "none") +
  theme(
    axis.title.y = element_markdown(face = "bold")
  ) +
  scale_color_manual(values = sample.cols)

# Fig. 2d
ggplot(depth.df %>% filter(!is.na(timepoints)), 
       aes(x = timepoints, y = log10(depth), color = timepoints)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 21, fill = NA, width = 0.1) +
  theme_minimal_hgrid() +
  labs(x = NULL, color = "Sample", y = "log~10~Depth") +
  guides(color = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme(
    axis.title.y = element_markdown(face = "bold"),
    strip.text = element_text(face = "bold")
  ) + scale_color_manual(values = tp.cols)

# Fig. 2e
ggplot(depth.df %>% filter(host_sex %in% c("M", "F")), 
       aes(x = host_sex, y = log10(depth), color = host_sex)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 21, fill = NA, width = 0.1) +
  theme_minimal_hgrid() +
  labs(x = NULL, color = "Sample", y = "log~10~Depth") +
  scale_color_manual(name = "Sex", values = sex.col) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  facet_wrap(vars(timepoints), nrow = 1) +
  guides(color = guide_legend(title.position="top", nrow = 1)) +
  theme(
    axis.title.y = element_markdown(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold")
  ) +
  stat_pwc(
    method = "wilcox_test",
    label = "{p.signif}"
  )

# Fig. 2f-h: Alpha rarefaction plot ----------------------------------------------------------------

# Processing metadata
md <- data.frame(sample_data(ps))
colnames(md)[1] <- "sample_id"
rpl.text <- c("Large_intestine" = "Li", "Rumen" = "Ru", "Small_intestine" = "Si")
md <- md %>% mutate(sample_origin = str_replace_all(sample_origin, rpl.text))

sex.lev <- c("M", "F")
day.level <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9")
sample.level <- c("Feces", "Li", "Si", "Ru",  "Milk", "Diet")
sex.col <- c(M = "blue", `F` = "red")
tp.cols <- c("#000000", "#004000", "#008000", "#00C000", "#01FE00", "#41BE00", 
             "#817E00", "#C13E00", "#FF6E00")
names(tp.cols) <- day.level
sample.cols <- c(Feces = "#307F00", Li = "#FFBF00", 
                 Si = "darkorange", Ru = "#FF0000", Milk = "#0F52BA", Diet = "violet")

# Alpha rarefaction is calculate from Qiime2 platform
alpha.rarefaction <- read_delim("Shannon_alpha_rarefaction.csv") # Load alpha rarefaction
colnames(alpha.rarefaction)[1] <- "sample_id"

# Processing data
mdf <- alpha.rarefaction %>% pivot_longer(!sample_id) %>% filter(!is.na(value))
mdf <- mdf %>% mutate(depth = str_extract(name, "^depth-([0-9]+)_", group = 1)) %>% mutate(depth = as.numeric(depth))
mdf <- mdf %>% group_by(sample_id, depth) %>% summarise(value = mean(value))
mdf <- mdf %>% left_join(md, by = "sample_id")
mdf$host_sex <- factor(mdf$host_sex, levels = sex.lev)
mdf$timepoints <- factor(mdf$timepoints, levels = day.level)
mdf$sample_origin <- factor(mdf$sample_origin, levels = sample.level)

# Fig. 2f
mdf %>% 
  group_by(sample_origin, depth) %>% summarise(mean = mean(value), 
                                               sd = sd(value), 
                                               n = n(),
                                               sem = sd/sqrt(n),
                                               t_sem_95 = qt(0.975, df = n-1)*sem) %>% ungroup %>% 
  ggplot(aes(x = depth, y = mean)) +
  geom_ribbon(aes(ymin = mean-t_sem_95, ymax = mean+t_sem_95, fill = sample_origin), alpha = 0.2) +
  geom_line(aes(group = sample_origin, color = sample_origin)) +
  scale_color_manual(name = "Sample type", values = sample.cols) +
  scale_fill_manual(name = "Sample type", values = sample.cols) +
  theme_cowplot() +
  labs(x = "Depth", y = "Shannon") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(
    legend.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

# Fig. 2g
mdf %>% filter(!is.na(timepoints)) %>% 
  group_by(timepoints, depth) %>% summarise(mean = mean(value), 
                                           sd = sd(value), 
                                           n = n(),
                                           sem = sd/sqrt(n),
                                           t_sem_95 = qt(0.975, df = n-1)*sem) %>% ungroup %>% 
  ggplot(aes(x = depth, y = mean)) +
  geom_ribbon(aes(ymin = mean-t_sem_95, ymax = mean+t_sem_95, fill = timepoints), alpha = 0.2) +
  geom_line(aes(group = timepoints, color = timepoints)) +
  scale_color_manual(name = "Timepoints", values = tp.cols) +
  scale_fill_manual(name = "Timepoints", values = tp.cols) +
  theme_cowplot() +
  labs(x = "Depth", y = "Shannon") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(
    legend.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

# Fig. 2h
mdf %>% filter(host_sex %in% c("M", "F")) %>% 
  group_by(host_sex, depth) %>% summarise(mean = mean(value), 
                                          sd = sd(value), 
                                          n = n(),
                                          sem = sd/sqrt(n),
                                          t_sem_95 = qt(0.975, df = n-1)*sem) %>% 
  ungroup() %>% 
  ggplot(aes(x = depth, y = mean)) +
  geom_ribbon(aes(ymin = mean-t_sem_95, ymax = mean+t_sem_95, fill = host_sex), alpha = 0.2) +
  geom_line(aes(group = host_sex, color = host_sex)) +
  scale_color_manual(values = sex.col) +
  scale_fill_manual(values = sex.col) +
  theme_cowplot() +
  labs(x = "Depth", y = "Shannon") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(
    legend.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )


# Fig. 3 ------------------------------------------------------------------

library(vegan)
library(rbiom)
library(rstatix)

beta.ps <- ps %>% subset_samples(!sample_origin %in% c("Diet", "Milk"))
rpl.text <- c("Large_intestine" = "Li", "Rumen" = "Ru", "Small_intestine" = "Si")
sample_data(beta.ps)$sample_origin <- str_replace_all(sample_data(beta.ps)$sample_origin, rpl.text)

tp.suffix <- str_replace(sample_data(beta.ps)$sample_origin, "Feces", "")
new.col <- paste(sample_data(beta.ps)$timepoints, tp.suffix, sep = "_")
new.col <- str_replace(new.col, "_$", "")
sample_data(beta.ps)$tp_so <- new.col

lev <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9_Li", "T9_Si", "T9_Ru")
cols <- c("#000000", "#004000", "#008000", "#00C000", "#01FE00", "#41BE00", 
          "#817E00", "#C13E00", "#FFBF00", "darkorange", "#FF0000")
names(cols) <- lev
sp <- c(rep(16, 8), 17, 15, 18)
names(sp) <- lev

min.depth <- floor(min(sample_sums(beta.ps)))
rf.ps <-  rarefy_even_depth(beta.ps, min.depth, replace = FALSE, rngseed = 123)
beta.md <- data.frame(rf.ps@sam_data)
asv.table <- data.frame(rf.ps@otu_table, check.names = FALSE)
dist <- unifrac(as.matrix(asv.table), weighted = TRUE, tree = phy_tree(rf.ps))
pcoa <- cmdscale(dist, k = nrow(as.matrix(dist))-1, eig = TRUE)

# PERMANOVA
sex.adonis <- adonis2(dist ~ host_sex, beta.md, permutations = 999)
ci.adonis <- adonis2(dist ~ calf_id, beta.md, permutations = 999)
tp.adonis <- adonis2(dist ~ timepoints, beta.md, permutations = 999)

ad.mat <- matrix(nrow = 3, ncol = 2)
ad.mat[,1] <- c(sex.adonis$R2[1], ci.adonis$R2[1], tp.adonis$R2[1])
ad.mat[,2] <-  c(sex.adonis$`Pr(>F)`[1], ci.adonis$`Pr(>F)`[1], tp.adonis$`Pr(>F)`[1])
ad.df <- data.frame(ad.mat)
names(ad.df) <- c("R2", "p")
ad.df <- add_significance(ad.df)
ad.df$var.name <- c("Sex", "Cattle ID", "Timepoints")

## Fig. 3c -----------------------------------------------------------------
ggplot(ad.df, aes(x = R2, y = fct_reorder(var.name, R2))) +
  geom_bar(stat = "identity", fill = "black") +
  labs(x = "R^2^", y = NULL) +
  geom_text(aes(x = max(R2)*0.1+R2, label = p.signif, y = as.numeric(fct_reorder(var.name, R2))), 
            fontface = "bold", size = 6) +
  theme_minimal_vgrid() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
    axis.title.x = element_markdown()
  )

# Draw PCoA plot
axis.mapping <- paste0("Axis", c(1, 2))
ordinations <- data.frame(pcoa$points, check.names = FALSE)
ordinations <- ordinations[,c(1, 2)]
colnames(ordinations) <- paste0("Axis", seq(1, ncol(ordinations)))
ordinations <- ordinations %>% rownames_to_column("sample_name")
ordi.df <- left_join(ordinations, beta.md, by = "sample_name")

# Make target column to factor
ordi.df$tp_so <- factor(ordi.df$tp_so, levels = lev)

# Calculate relative eigen values
eig <- pcoa$eig
rel.eig <- eig / sum(eig)

# Axis parameters
axis.labels <- round(rel.eig[c(1, 2)]*100, 2)
axis.labels <- as.character(glue("PCo{c(1, 2)} ({axis.labels}%)"))
segment.mapping <- paste0("end_", axis.mapping)

## Fig. 3a -----------------------------------------------------------------

ggplot(ordi.df, aes(x = !!sym(axis.mapping[1]), y = !!sym(axis.mapping[2]))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = .5) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = .5) +
  stat_ellipse(
    aes(group = tp_so, color = tp_so), 
    lwd = 0.5, linetype = 2, show.legend = FALSE, segments = 101
  ) +
  geom_point(
    aes(color = tp_so, shape = tp_so),
    size = 2.5) + theme_bw() + 
  labs( 
    x = axis.labels[1], y = axis.labels[2]
  ) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.text = element_markdown(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal"
  ) +  
  scale_color_manual(name = "Timepoints", values = cols) +
  scale_shape_manual(name = "Timepoints", values = sp)


library(ggrepel)

axis <- c(1, 2)
sp <- c(16, 17, 15, 18)
names(sp) <- c("Feces", "Li", "Si", "Ru")
axis.mapping <- paste0("Axis", axis)

ordinations <- data.frame(pcoa$points, check.names = FALSE)
ordinations <- ordinations[axis]
colnames(ordinations) <- paste0("Axis", seq(1, ncol(ordinations)))
ordinations <- ordinations %>% rownames_to_column("sample_name")
ordi.df <- left_join(ordinations, beta.md, by = "sample_name")

cent <- ordi.df %>% group_by(tp_so) %>% 
  summarise(across(all_of(axis.mapping), ~mean(.x))) %>% 
  rename_with(~paste0("end_", .x), starts_with("Axis"))

ordi.df <- ordi.df %>% left_join(cent, by = "tp_so")

# Calculate relative eigen values
eig <- input.ordi$eig
rel.eig <- eig / sum(eig)

# Axis parameters
axis.labels <- round(rel.eig[axis]*100, 2)
axis.labels <- as.character(glue("PCo{axis} ({axis.labels}%)"))
segment.mapping <- paste0("end_", axis.mapping)
ordi.df <- ordi.df %>% arrange(tp_so)
ordi.df.mean.total <- ordi.df %>% 
  group_by(timepoints) %>% 
  summarise(!!sym(axis.mapping[1]) := mean(!!sym(axis.mapping[1])),
            !!sym(axis.mapping[2]) := mean(!!sym(axis.mapping[2]))) %>% 
  mutate(group = 1)

ordi.df.mean.sex <- ordi.df %>% 
  group_by(timepoints, host_sex) %>% 
  summarise(!!sym(axis.mapping[1]) := mean(!!sym(axis.mapping[1])),
            !!sym(axis.mapping[2]) := mean(!!sym(axis.mapping[2])))

## Fig. 3d-------------------------------------------------------------

ggplot(ordi.df, aes(x = !!sym(axis.mapping[1]), y = !!sym(axis.mapping[2]))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = .5) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = .5) + 
  stat_ellipse(data = ordi.df %>% filter(sample_origin != "Feces"),
                      aes(group = sample_origin), 
                      lwd = 0.5, linetype = 2, show.legend = FALSE, segments = 101) +
  geom_point(aes(shape = sample_origin), size = 2.5, alpha = 0.2, color = "black") +
  scale_shape_manual(name = "Sample type", values = sp) +
  guides(shape = guide_legend(override.aes = list(alpha = 1))) + 
  geom_path(data = ordi.df.mean.total, aes(group = group), linewidth = 1, 
                   color = "purple", 
                   arrow = arrow(length = unit(0.1, "inches"))) +
  geom_point(data = ordi.df.mean.total, shape = 21, size = 2.5, color = "purple", fill = "white") +
  geom_text_repel(data = ordi.df.mean.total, 
                         aes(label = timepoints), size = 5, color = "black") + theme_bw() + 
  labs(x = axis.labels[1], y = axis.labels[2]) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )


## Fig. 3e -----------------------------------------------------------------

ggplot(ordi.df, aes(x = !!sym(axis.mapping[1]), y = !!sym(axis.mapping[2]))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = .5) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = .5) + 
  stat_ellipse(data = ordi.df %>% filter(sample_origin != "Feces"),
                      aes(group = interaction(sample_origin, host_sex), color = host_sex), 
                      lwd = 0.5, linetype = 2, show.legend = FALSE, segments = 101) +
  geom_point(aes(color = host_sex, shape = sample_origin),size = 2.5, alpha = 0.2) +
  scale_shape_manual(name = "Sample type", values = sp) + 
  geom_path(data = ordi.df.mean.sex, aes(color = host_sex, group = host_sex), linewidth = 1, 
                   arrow = arrow(length = unit(0.1, "inches")),
                   key_glyph = "rect") +
  geom_point(data = ordi.df.mean.sex, aes(color = host_sex), size = 2.5, shape = 21, fill = "white",
             show.legend = FALSE) + 
  geom_text_repel(data = ordi.df.mean.sex, aes(label = timepoints, color = host_sex),
                         box.padding = 0.1, size = 3.5, show.legend = FALSE) + theme_bw() + 
  labs(x = axis.labels[1], y = axis.labels[2]) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.text = element_markdown(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.box = "vertical",
    legend.key.size = unit(5, "mm")
  ) + scale_color_manual(values = c(M = "blue", `F` = "red")) +
  guides(color = guide_legend(order=1), 
         shape = guide_legend(order=2, override.aes = list(alpha = 1)))

## Fig. 3b -----------------------------------------------------------------

beta.ps <- ps %>% subset_samples(sample_origin %in% c("Large_intestine", "Small_intestine", "Rumen"))
rpl.text <- c("Large_intestine" = "Li", "Rumen" = "Ru", "Small_intestine" = "Si")
sample_data(beta.ps)$sample_origin <- str_replace_all(sample_data(beta.ps)$sample_origin, rpl.text)

cols <- c(Li="#FFBF00", Si="darkorange", Ru="#FF0000")
sp <- c(Li=17, Si=15, Ru=18)

min.depth <- floor(min(sample_sums(beta.ps)))
rf.ps <-  rarefy_even_depth(beta.ps, min.depth, replace = FALSE, rngseed = 123)
beta.md <- data.frame(rf.ps@sam_data)
asv.table <- data.frame(rf.ps@otu_table, check.names = FALSE)
dist <- unifrac(as.matrix(asv.table), weighted = TRUE, tree = phy_tree(rf.ps))
pcoa <- cmdscale(dist, k = nrow(as.matrix(dist))-1, eig = TRUE)

# Draw PCoA plot
axis.mapping <- paste0("Axis", c(1, 2))
ordinations <- data.frame(pcoa$points, check.names = FALSE)
ordinations <- ordinations[,c(1, 2)]
colnames(ordinations) <- paste0("Axis", seq(1, ncol(ordinations)))
ordinations <- ordinations %>% rownames_to_column("sample_name")
ordi.df <- left_join(ordinations, beta.md, by = "sample_name")

# Calculate relative eigen values
eig <- pcoa$eig
rel.eig <- eig / sum(eig)

# Axis parameters
axis.labels <- round(rel.eig[c(1, 2)]*100, 2)
axis.labels <- as.character(glue("PCo{c(1, 2)} ({axis.labels}%)"))
segment.mapping <- paste0("end_", axis.mapping)

ggplot(ordi.df, aes(x = !!sym(axis.mapping[1]), y = !!sym(axis.mapping[2]))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = .5) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = .5) +
  stat_ellipse(
    aes(group = sample_origin, color = sample_origin), 
    lwd = 0.5, linetype = 2, show.legend = FALSE, segments = 101
  ) +
  geom_point(
    aes(color = sample_origin, shape = sample_origin),
    size = 2.5) + theme_bw() + 
  labs( 
    x = axis.labels[1], y = axis.labels[2]
  ) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank()
  ) +  
  scale_color_manual(name = "Timepoints", values = cols) +
  scale_shape_manual(name = "Timepoints", values = sp)

## Fig. 3f ----------------------------------------------------

day.level <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9_Li", "T9_Si", "T9_Ru", "T9")
cols <- c(M = "blue", `F` = "red")
sample.level <- c("Feces", "Large_intestine", "Small_intestine", "Rumen",  "Milk", "Diet")

# Day_order_sample
for (cur.day in day.level){
  
  if (cur.day == "T9"){
    cur.ps <- beta.ps %>% subset_samples(timepoints == cur.day)
  } else {
    cur.ps <- beta.ps %>% subset_samples(tp_so == cur.day)  
  }
  
  min.depth <- floor(min(sample_sums(cur.ps)))
  rf.ps <-  rarefy_even_depth(cur.ps, min.depth, replace = FALSE, rngseed = 123)
  beta.md <- data.frame(rf.ps@sam_data)
  asv.table <- data.frame(rf.ps@otu_table, check.names = FALSE)
  dist <- unifrac(as.matrix(asv.table), weighted = TRUE, tree = phy_tree(rf.ps))
  pcoa <- cmdscale(dist, k = nrow(as.matrix(dist))-1, eig = TRUE)
  
  # PERMANOVA
  res <- adonis2(dist ~ host_sex, beta.md, permutations = 999)

  axis <- c(1, 2)
  axis.mapping <- paste0("Axis", axis)
  ordinations <- data.frame(pcoa$points, check.names = FALSE)
  ordinations <- ordinations[axis]
  colnames(ordinations) <- paste0("Axis", seq(1, ncol(ordinations)))
  ordinations <- ordinations %>% rownames_to_column("sample_name")
  ordi.df <- left_join(ordinations, beta.md, by = "sample_name")

  # Calculate relative eigen values
  eig <- pcoa$eig
  rel.eig <- eig / sum(eig)
  
  # Axis parameters
  axis.labels <- round(rel.eig[axis]*100, 2)
  axis.labels <- as.character(glue("PCo{axis} ({axis.labels}%)"))
  segment.mapping <- paste0("end_", axis.mapping)
  
  ggplot(ordi.df, aes(x = !!sym(axis.mapping[1]), y = !!sym(axis.mapping[2]))) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = .5) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = .5) +
    stat_ellipse(
      aes(group = host_sex, color = host_sex), 
      lwd = 0.5, linetype = 2, show.legend = FALSE, segments = 101
    ) +
    geom_point(
      aes(color = host_sex),
      size = 2.5) +
    annotate(geom = "text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, label = cur.day) + theme_bw() + 
    labs(x = axis.labels[1], y = axis.labels[2]) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      legend.text = element_markdown(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.box = "horizontal"
    ) + scale_color_manual(name = "Sex", values = cols)
}




# Fig. 4 ------------------------------------------------------------------

alpha.ps <- ps %>% subset_samples(!sample_origin %in% c("Diet", "Milk"))
min.depth <- floor(min(sample_sums(alpha.ps)))
rf.ps <-  rarefy_even_depth(alpha.ps, min.depth, replace = FALSE, rngseed = 123)

alphaDiversity <- estimate_richness(rf.ps, measures = "Shannon") 
alpha.df <- left_join(alphaDiversity %>% rownames_to_column("merge.id"), 
                      data.frame(rf.ps@sam_data) %>% rownames_to_column("merge.id"), 
                      by = "merge.id")
alpha.df <- alpha.df %>% column_to_rownames("merge.id")
alpha.mdf <- alpha.df %>% pivot_longer(Shannon, names_to = "Alpha_measures", values_to = "Alpha_div_values")

rpl.text <- c("Large_intestine" = "Li", "Rumen" = "Ru", "Small_intestine" = "Si")
alpha.mdf <- alpha.mdf %>% mutate(sample_origin = str_replace_all(sample_origin, rpl.text))

tp.suffix <- str_replace(alpha.mdf$sample_origin, "Feces", "")
new.col <- paste(alpha.mdf$timepoints, tp.suffix, sep = "_")
new.col <- str_replace(new.col, "_$", "")
alpha.mdf$tp_so <- new.col

alpha.mdf$tp_so <- factor(alpha.mdf$tp_so, levels = c("T1","T2","T3","T4","T5","T6","T7","T8","T9_Li","T9_Si","T9_Ru"))
alpha.mdf$timepoints <- factor(alpha.mdf$timepoints, levels = c("T1","T2","T3","T4","T5","T6","T7","T8","T9"))
alpha.mdf$host_sex <- factor(alpha.mdf$host_sex, levels = c("M", "F"))

alpha.mdf.t9 <- alpha.mdf %>% filter(timepoints == "T9") %>% 
  group_by(calf_id, timepoints, host_sex) %>% 
  summarise(Alpha_div_values = mean(Alpha_div_values))

alpha.mdf.day <- bind_rows(alpha.mdf %>% filter(timepoints != "T9"), alpha.mdf.t9)


## Fig. 4a -----------------------------------------------------------------

# Fig. 4a left panel
ggplot(alpha.mdf.day, aes(x = timepoints, y = Alpha_div_values)) +
  geom_point(color = "lightgray", size = 0.75) +
  geom_line(aes(group = calf_id), color = "lightgray") +
  stat_summary(geom = "point", fun = "mean") +
  stat_summary(geom = "line", mapping = aes(group = 1), fun = "mean") +
  theme_cowplot() +
  labs(x = "Timepoints", y = "Shannon") +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Fig. 4a right panel
alpha.mdf %>% filter(sample_origin != "Feces") %>% 
  ggplot(aes(x = tp_so, y = Alpha_div_values)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.05) +
  theme_cowplot() +
  labs(x = "T9_samples", y = NULL) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_continuous(limits = y.limits)


## Fig. 4b -----------------------------------------------------------------

# Fig. 4b left panel
ggplot(alpha.mdf.day, aes(x = timepoints, y = Alpha_div_values)) +
  geom_line(aes(group = calf_id, color = host_sex), alpha = 0.2, linewidth = 0.4, key_glyph = "rect") +
  geom_point(aes(color = host_sex), alpha = 0.2, size = 0.75, key_glyph = "rect") +
  stat_summary(geom = "line", mapping = aes(group = host_sex, color = host_sex), fun = "mean", 
               linewidth = 0.75, key_glyph = "rect") +
  stat_summary(geom = "point", mapping = aes(color = host_sex), fun = "mean",
               size = 1.5, shape = 21, fill = "white", show.legend = FALSE) +
  theme_cowplot() +
  labs(x = "Timepoints", y = "Shannon") +
  scale_color_manual(name = "Sex", values = c(M = "blue", `F` = "red")) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Fig. 4b right panel
alpha.mdf %>% filter(sample_origin != "Feces") %>% 
  ggplot(aes(x = tp_so, y = Alpha_div_values, color = host_sex)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.05)) +
  theme_cowplot() +
  labs(x = "T9_samples", y = NULL) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = y.limits) +
  scale_color_manual(values = c(M = "blue", `F` = "red"))

# Wilcoxon rank-sum test with Bonferroni–Holm correction
alpha.mdf %>% 
  group_by(tp_so) %>% 
  wilcox_test(Alpha_div_values ~ host_sex) %>% 
  adjust_pvalue(method = "holm") %>% 
  add_significance()


# Fig. 5 ------------------------------------------------------------------

library(ggalluvial)
library(microbiome)
library(RColorBrewer)

input.ps <- ps
lev1 <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9_Li", "T9_Si", "T9_Ru", "Diet", "Milk")
cols <- c("#000000", "#004000", "#008000", "#00C000", "#01FE00", "#41BE00", 
          "#817E00", "#C13E00", "#FFBF00", "darkorange", "#FF0000", "#0F52BA", "violet")
names(cols) <- lev1

other.rpl <- "Others(≤5%RA)"

# Choose one of the taxa level below.
# Use Phylum level for Fig. 5a, b and Genus level for Fig. 5c, d.
tax_level <- c("Phylum", "Genus")

agg.ps <- input.ps %>% microbiome::transform("compositional") %>% 
  aggregate_rare(level = tax_level, detection = 0.05, prevalence = 0)
mdf <- psmelt2(agg.ps)
mdf <- mdf %>% mutate(FeatureID = str_replace(FeatureID, "Other", other.rpl))
rpl.text <- c("Large_intestine" = "Li", "Rumen" = "Ru", "Small_intestine" = "Si")
mdf <- mdf %>% mutate(sample_origin = str_replace_all(sample_origin, rpl.text))

tp.suffix <- str_replace(mdf$sample_origin, "Feces", "")
new.col <- paste(mdf$timepoints, tp.suffix, sep = "_")
new.col <- str_replace(new.col, "_$", "")
mdf$tp_so <- new.col

mdf$tp_so <- factor(mdf$tp_so, levels = lev1)
mdf.total <- mdf %>% group_by(FeatureID, tp_so) %>% summarise(value = mean(value)) %>% ungroup

mdf.sex <- mdf %>% group_by(FeatureID, tp_so, host_sex) %>% summarise(value = mean(value)) %>% ungroup
mdf.sex <- mdf.sex %>% filter(host_sex %in% c("M", "F"))
mdf.sex$host_sex <- factor(mdf.sex$host_sex, levels = c("M", "F"))

# Set taxa level
abn_level <- mdf.total %>% group_by(FeatureID) %>% summarise(value = sum(value)) %>% 
  arrange(desc(value)) %>% pull(FeatureID) %>% unique
other.idx <- which(abn_level == other.rpl)
if (!is_empty(other.idx)) abn_level <- c(abn_level[-other.idx], abn_level[other.idx])

mdf.total$FeatureID <- factor(mdf.total$FeatureID, levels = abn_level)
mdf.sex$FeatureID <- factor(mdf.sex$FeatureID, levels = abn_level)

# Set taxa colors
x <- unique(mdf.total$FeatureID)
n <- length(unique(x))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
is.pastel <- str_detect(rownames(qual_col_pals), "Pastel")
qual_col_pals = bind_rows(qual_col_pals[!is.pastel,], qual_col_pals[is.pastel,])
is.set <- str_detect(rownames(qual_col_pals), "Set")
qual_col_pals <- bind_rows(qual_col_pals[is.set,], qual_col_pals[!is.set,])
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
mycolors <- rep_len(col_vector, length.out = n)
names(mycolors) <- unique(x)
mycolors[other.rpl] <- "black"


## Fig. 5a, c --------------------------------------------------------------

p <- ggplot(mdf.total, aes(x = tp_so, y = value, 
                           fill = FeatureID)) + 
  geom_stratum(aes(stratum = FeatureID), color = "black", width = 3/5, size = 0.15) +
  geom_flow(aes(alluvium = FeatureID), curve_type = "linear", width = 3/5) +
  labs(y = "Relative abundance", x = NULL) +
  theme_minimal_hgrid() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(4, 'mm'),
    legend.spacing.y = unit(0.5, "mm")
  ) + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) +
  guides(fill = guide_legend(title = tax_level, byrow = TRUE))

## Fig. 5b, d --------------------------------------------------------------
p <- ggplot(mdf.sex, aes(x = host_sex, y = value, fill = FeatureID)) + 
  facet_wrap(vars(tp_so), scales = "free_x", nrow = 1) +
  geom_stratum(aes(stratum = FeatureID), color = "black", width = 3/5, size = 0.15) +
  geom_flow(aes(alluvium = FeatureID), curve_type = "linear", width = 3/5) +
  labs(y = "Relative abundance", x = NULL) +
  theme_minimal_hgrid() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 11),
    strip.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(4, 'mm'),
    legend.spacing.y = unit(0.5, "mm")
  ) + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors)

