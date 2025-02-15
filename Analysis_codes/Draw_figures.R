#################################################################################################
########## Assembly and maturation of calf gut microbiome from neonate to post-puberty ##########
#################################################################################################
# Code written by Jae-Yun Lee
# Revised at 2024-07-02

library(tidygraph)
library(igraph)
library(microbiome)
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
tax.df <- tax.df %>% column_to_rownames("ASV")
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


# Figure 3: Top 20 prevalent taxonomy --------------------------------------------------------------------

n.samples <- nsamples(ps)
prv.cut <- n.samples*0.1
agg.ps <- aggregate_taxa(ps, "Genus")

filt.ps <- transformSampleCounts(agg.ps, function(.x) ifelse(.x > 0, 1, 0))
taxa.prev <- sort(taxa_sums(filt.ps), decreasing = T)
taxa.prev <- taxa.prev / n.samples * 100
target.tax <- head(taxa.prev, 20)

tax.df <- data.frame(filt.ps@tax_table)
tax.df <- tax.df[names(target.tax),]
rownames(tax.df) <- NULL

nodes <- tax.df %>%
  pivot_longer(cols = Kingdom:Genus, names_to = "level", values_to = "taxa") %>%
  group_by(unique) %>%
  mutate(Parent = lag(taxa, default = NA)) %>%
  mutate(Phylum = taxa[level == "Phylum"]) %>% ungroup %>% 
  ungroup()

graph.df <- nodes %>% 
  filter(!is.na(Parent)) %>%
  select(Parent, taxa, Phylum, level) %>%
  mutate(Phylum = ifelse(Parent == "Bacteria", NA, Phylum)) %>% 
  distinct()

graph <- graph_from_data_frame(graph.df, directed = TRUE)

V(graph)$level <- sapply(V(graph)$name, function(x) {
  nodes %>% filter(taxa == x) %>% pull(level) %>% unique
}) %>% unlist %>% unname
phylum.vec <- nodes$Phylum
names(phylum.vec) <- nodes$taxa
phylum.vec[names(phylum.vec) == "Bacteria"] <- NA
V(graph)$Phylum <- phylum.vec[V(graph)$name]

length.info <- c(seq(1, 5), 6.5) 
names(length.info) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
V(graph)$height <- length.info[V(graph)$level]

ggraph(graph, layout = "dendrogram", height = height) +  
  geom_edge_elbow(aes(color = Phylum), edge_width = 1) +  
  geom_node_point(aes(color = Phylum), shape = 21, size = 6, fill = "white", stroke = 1) + 
  geom_node_label(aes(label = name), vjust = 1.5, size = 4, alpha = 0.5, label.size = NA) +
  scale_x_discrete(expand = expansion(mult = 0.5)) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_edge_color_manual(values = c(Bacteroidota = "darkgreen", Firmicutes = "darkblue", Proteobacteria = "darkorange"),
                          na.value = "black") +
  scale_color_manual(values = c(Bacteroidota = "darkgreen", Firmicutes = "darkblue", Proteobacteria = "darkorange"),
                     na.value = "black") +
  theme_void() +
  theme(legend.position = "none")
