#-----------------------------------------------------
# Title: Differential Abundance
# Date: 11/08/2023
# Author: Arif Ayman
# Description: Differential Abundance analysis with Aldex2
#-----------------------------------------------------


# Generate dataframe with OTUs and metadata
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox$Group <- phyloseq::sample_data(ps_clr)$Group

# Define functions to pass to map
wilcox_model <- function(df) {
  wilcox.test(abund ~ Group, data = df)
}
wilcox_pval <- function(df) {
  wilcox.test(abund ~ Group, data = df)$p.value
}

# Create nested data frames by OTU and loop over each using map 
wilcox_results <- ps_wilcox %>%
  gather(key = OTU, value = abund, -Group) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))  

# Show results
head(wilcox_results)

head(wilcox_results$data[[1]])

wilcox_results$wilcox_test[[1]]

wilcox_results$p_value[[1]]


# Unnesting
wilcox_results <- wilcox_results %>%
  dplyr::select(OTU, p_value) %>%
  unnest(cols = c(p_value))

#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps_clr))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
# Modify OTU values in wilcox_results
wilcox_results$OTU <- str_replace(wilcox_results$OTU, "OTU\\.", "OTU ")

# Join with taxa_info and include taxonomic information
wilcox_results <- wilcox_results %>%
  left_join(taxa_info, by = "OTU") %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(OTU, p_value, BH_FDR, Kingdom:Species, everything())

print.data.frame(wilcox_results)



# Run ALDEx2 with Wilcoxon test
aldex2_da <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$Group, test = "t", effect = TRUE, denom = "iqlr")

ALDEx2::aldex.plot(aldex2_da, type = "MW", test = "wilcox", called.cex = 1, cutoff = 0.05)


# Clean up presentation
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.82 ) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# Join with taxa_info
sig_aldex2 <- left_join(sig_aldex2, taxa_info, by = "OTU")

sig_aldex2


# Identify species more abundant in Definite_IE group
definite_ie_species <- sig_aldex2 %>%
  filter(diff.btw > 0) %>%
  pull(Species)

# Identify species more abundant in IE_High_Risk group
high_risk_species <- sig_aldex2 %>%
  filter(diff.btw < 0) %>%
  pull(Species)

# Print the results
if (length(definite_ie_species) > 0) {
  cat(paste(definite_ie_species, "species is more abundant in Definite_IE.", "\n"))
} else {
  cat("No species are more abundant in Definite_IE.", "\n")
}

if (length(high_risk_species) > 0) {
  cat(paste(high_risk_species, "species is more abundant in IE_High_Risk.", "\n"))
} else {
  cat("No species are more abundant in IE_High_Risk.", "\n")
}
