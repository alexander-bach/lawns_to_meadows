###############################################################

#From lawns to meadows: Spiders (Arachnida: Araneae) as indicators to 
#measure urban grassland restoration success

#Bach et al. 2024, Urban Ecosystems


###############################################################

# R code to reproduce results of the MRPP and the chao1 analysis.

# Author: Alexander Bach
# Last update: 01.09.2024

###############################################################

library(here)
library(data.table)
library(dplyr)
library(vegan)

#### load and prepare necessary data ####
data <- read.csv2(here("spider_data.csv"), sep = "\t", encoding = "UTF-8") 

sites <- read.csv2(here("sites.csv"), sep = "\t", encoding = "UTF-8") 
sites <- sites[ order(sites$LocalityDescription), ]

# sum up sampledata from different Collection timepoints
data <- group_by(data, LastIdentificationCache, LocalityDescription_2, CollectionTimeSpan, Value)
data <- summarise(data, NumberOfUnits = sum(NumberOfUnits))

# create crosstable for analyses
data_cr <- setDT(data) %>%
  dcast(LocalityDescription_2~LastIdentificationCache, value.var = "NumberOfUnits", sum)

data_cr <- as.data.frame(data_cr)
row.names(data_cr) <- data_cr$LocalityDescription_2
data_cr$LocalityDescription_2 <- NULL

data_cr <- data_cr[ order(row.names(data_cr)), ]

# merge data for analyses
merged_data <- merge(data_cr, sites[, c("LocalityDescription", "Naturalness")], 
                     by.x = "row.names", by.y = "LocalityDescription")
rownames(merged_data) <- merged_data$Row.names
merged_data$Row.names <- NULL
species_data <- merged_data[, !colnames(merged_data) %in% c("Naturalness")]
groups <- merged_data$Naturalness

##### calculate MRPP ####
pairwise_mrpp <- function(species_data, groups) {
  unique_groups <- unique(groups)
  results <- data.frame(Group1 = character(), Group2 = character(), p.value = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:(length(unique_groups) - 1)) {
    for (j in (i + 1):length(unique_groups)) {
      group1 <- unique_groups[i]
      group2 <- unique_groups[j]
      
      # Subset data for the pair of groups
      subset_data <- species_data[groups %in% c(group1, group2), ]
      subset_groups <- groups[groups %in% c(group1, group2)]
      
      # Perform MRPP
      mrpp_result <- mrpp(subset_data, subset_groups, permutations = 999)
      
      # Store the result
      results <- rbind(results, data.frame(Group1 = group1, Group2 = group2, p.value = mrpp_result$Pvalue))
    }
  }
  
  return(results)
}

# Perform pairwise MRPP
pairwise_results <- pairwise_mrpp(species_data, groups)

# Adjust p-values using Holm method
pairwise_results$adjusted_p.value <- p.adjust(pairwise_results$p.value, method = "holm")

# Print the pairwise results
print(pairwise_results)

#### calculate chao1 ####

# Aggregate species data by groups
species_by_group <- aggregate(. ~ groups, data = data.frame(species_data, groups), FUN = sum)

# Set row names as groups
rownames(species_by_group) <- species_by_group$groups
species_by_group$groups <- NULL

# Calculate Chao1 index, Chao1 SE, and observed species richness for each group
chao1_results <- apply(species_by_group, 1, function(x) {
  est <- estimateR(x)
  list(S.chao1 = est[["S.chao1"]], SE.chao1 = est[["se.chao1"]])
})

# Extract Chao1 index and SE into separate vectors
chao1_index <- sapply(chao1_results, `[[`, "S.chao1")
chao1_se <- sapply(chao1_results, `[[`, "SE.chao1")

# Calculate observed species richness
observed_species <- rowSums(species_by_group > 0)

# Combine Chao1 index, Chao1 SE, and observed species richness into a data frame
results_df <- data.frame(
  Group = names(chao1_results),
  Observed_Species = observed_species,
  Chao1_Index = chao1_index,
  Chao1_SE = chao1_se
)