#-----------------------------------------------------
# Title: Beta Diversity 
# Date: 11/08/2023
# Author: Arif Ayman
# Description: Beta Diversity and Ordination analysis
#-----------------------------------------------------

ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
as.matrix(bc_dist)[1:5, 1:5]

# Dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))

# Color code
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
colorCode <- c(Definite_IE = "red", `IE_High_Risk` = "blue")
labels_colors(ward) <- colorCode[meta$Group][order.dendrogram(ward)]

# Plotting the dendrogram with improved readability
plot(ward, main = "Bray-Curtis Dissimilarity of Samples",
     xlab = "", ylab = "", sub = "", cex = 0.8)




# CLR transformation
(ps_clr <- microbiome::transform(ps, "clr"))
# PCA
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
# Scree plot
phyloseq::plot_scree(ord_clr) +
  labs(x = "\nAxis", y = "Proportion of Variance\n", title = "Principal Component Analysis")
# Scale and plot
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps, ord_clr, type = "SampleID", color = "Group", title = "Proportion of Total Variation") +
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Group), linetype = 2)

# Distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean")

# ADONIS
vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$Group)


# Dispersion and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$Group)
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

# Distance to Centroid
boxplot(dispr, main = "Centroid Distance Comparison", xlab = "", col = c("Definite_IE" = "#CC001B", "IE_High_Risk" = "#7D92E7")) 

# Permutation test
permutest(dispr)


# Combine disease and control samples again
ps <- merge_phyloseq(disease_samples, control_samples)

ps_ra <- transform_sample_counts(ps, function(x) x/sum(x))

# Calculate distance matrices
distance_mat <- phyloseq::distance(ps_ra, method = "bray")

# Plotting using Principal Coordinate Analysis (PCoA)
pcoa <- ordinate(ps, method = "PCoA", distance = distance_mat)

pcoa.plot <- plot_ordination(ps, pcoa, color = "Group", shape = "Group") +
  geom_point(size = 4, aes(fill = Group)) + # Fill color based on Group, adjust size
  scale_color_brewer(palette = "Set1") + # Change color palette
  scale_shape_manual(values=c(21, 22, 23)) + # Change point shapes
  geom_text_repel(aes(label = SampleID)) + # Add sample labels
  theme_bw() + # Use black & white theme for a clean look
  labs(
    title = "PCoA - Bray-Curtis",
    color = "Group",
    shape = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14), # Make the title centered, bold, and larger
    legend.position = "bottom" # Position the legend at the bottom
  )

print(pcoa.plot)

# Run PERMANOVA test
permanova_result <- vegan::adonis2(distance_mat ~ sample_data(ps)$Group, permutations = 999)

# Print the results
print(permanova_result)

