#-----------------------------------------------------
# Title: Relative Abundance 
# Date: 11/08/2023
# Author: Arif Ayman
# Description: Relative abundance graphs 
#-----------------------------------------------------

# Removing low abundance taxa
ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)

# Set Group as factor
sample_data(ps)$Group <- factor(sample_data(ps)$Group, levels = c("Definite_IE","IE_High_Risk"))

# Counting samples in each Group
sample_data(ps) %>% dplyr::count(Group)

# Phyla Count
table(tax_table(ps)[, "Phylum"])

# Converting to relative abundance
ps_rel_abund = transform_sample_counts(ps, function(x) x / sum(x))

# Replace long labels with shorter aliases in tax table
tax_table(ps_rel_abund)[tax_table(ps_rel_abund)[, "Phylum"] == "Candidatus_Saccharibacteria", "Phylum"] <- "C. Saccharibacteria"
tax_table(ps_rel_abund)[tax_table(ps_rel_abund)[, "Phylum"] == "Candidatus_Absconditabacteria", "Phylum"] <- "C. Absconditabacteria"

plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  ggtitle("Relative Abundance of Phyla") +
  facet_wrap(~ Group, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"))




# Agglomerating to phylum and renaming
ps_phylum <- tax_glom(ps_rel_abund, "Phylum")
taxa_names(ps_phylum) <- tax_table(ps_phylum)[, "Phylum"]
taxa_names(ps_phylum)[taxa_names(ps_phylum) == "Candidatus_Saccharibacteria"] <- "C. Saccharibacteria"
taxa_names(ps_phylum)[taxa_names(ps_phylum) == "Candidatus_Absconditabacteria"] <- "C. Absconditabacteria"


# Melt and plot
psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = Group, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  ggtitle("Phylum-Level Abundance in IE High Risk and Definite IE") +
  facet_wrap(~ OTU, scales = "free") +
  theme(
    plot.title = element_text(size = 16), 
    strip.text = element_text(size = 12, face = "bold"), # Increase size and make facet titles bold
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) # Adjust the position and angle
