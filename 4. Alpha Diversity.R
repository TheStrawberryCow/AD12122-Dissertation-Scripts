#-----------------------------------------------------
# Title: Alpha Diversity 
# Date: 11/08/2023
# Author: Arif Ayman
# Description: Alpha Diverity analysis with multiple metrics
#-----------------------------------------------------

# Compute alpha diversity for Shannon index
alpha_div <- phyloseq::estimate_richness(ps, split = TRUE, measures = c("Shannon"))
alpha_div$SampleID <- rownames(alpha_div)
alpha_div$Group <- sample_data(ps)$Group[match(rownames(alpha_div), sample_names(ps))]

# Reorder the levels of SampleID based on Group
alpha_div$SampleID <- factor(alpha_div$SampleID, levels = alpha_div$SampleID[order(alpha_div$Group)])

# Prepare data for Shannon metric
shannon_div <- alpha_div[, c("SampleID", "Shannon", "Group")]

# Plot Shannon diversity
ggplot(shannon_div, aes(x = SampleID, y = Shannon, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) +
  labs(x = "Sample ID", y = "Alpha Diversity", title = "Shannon Alpha Diversity Analysis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) + 
  scale_fill_manual(values = c("Definite_IE" = "#CC001B", "IE_High_Risk" = "#7D92E7")) 



# Subset data
disease_samples <- subset_samples(ps, Group == "Definite_IE")
control_samples <- subset_samples(ps, Group == "IE_High_Risk")

# Disease
alpha_div_disease <- estimate_richness(na.omit(disease_samples), split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
alpha_div_disease$SampleID <- rownames(alpha_div_disease)
alpha_div_disease$Group <- "Definite_IE"

# High_Risk
alpha_div_control <- estimate_richness(na.omit(control_samples), split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
alpha_div_control$SampleID <- rownames(alpha_div_control)
alpha_div_control$Group <- "High_Risk"

# Combine data
alpha_div <- rbind(alpha_div_disease, alpha_div_control)

# Print the alpha diversity values
print(alpha_div)
mean_alpha_div <- aggregate(cbind(Observed, Shannon, Simpson) ~ Group, alpha_div, mean)
print(mean_alpha_div)


# Melt the data for ggplot
mean_alpha_div_melt <- melt(mean_alpha_div, id.vars = "Group")
# Create subsets of the data for each metric
mean_alpha_div_melt_observed <- mean_alpha_div_melt[mean_alpha_div_melt$variable == "Observed", ]
mean_alpha_div_melt_shannon <- mean_alpha_div_melt[mean_alpha_div_melt$variable == "Shannon", ]
mean_alpha_div_melt_simpson <- mean_alpha_div_melt[mean_alpha_div_melt$variable == "Simpson", ]

# Create function for the plot
plot_alpha_div <- function(data, metric) {
  ggplot(data, aes(x = variable, y = value, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    labs(x = "Alpha Diversity", y = paste("Mean", metric, "Index Value"), title = paste("Mean", metric, "Alpha Diversity Comparison")) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    scale_fill_manual(values = c("#CC001B", "#7D92E7")) +
    guides(fill = guide_legend(title = "Group"))
}

# Plot each metric separately
plot_alpha_div(mean_alpha_div_melt_observed, "Observed")
plot_alpha_div(mean_alpha_div_melt_shannon, "Shannon")
plot_alpha_div(mean_alpha_div_melt_simpson, "Simpson")



# Perform Wilcoxon test for Observed diversity
wilcox.observed <- wilcox.test(alpha_div$Observed ~ alpha_div$Group)


# Print the test results
cat("Wilcoxon test for Observed diversity:\n")
print(wilcox.observed)


#Richness with breakaway
ba_adiv <- breakaway::breakaway(ps)
ba_adiv[1]
plot(ba_adiv, ps, color = "Group", title = "Estimate richness of Samples")   
summary(ba_adiv) %>%
  add_column("SampleNames" = ps %>% otu_table %>% sample_names)  
bt <- breakaway::betta(summary(ba_adiv)$estimate,
                       summary(ba_adiv)$error,
                       make_design_matrix(ps, "Group"))
bt$table                                 


