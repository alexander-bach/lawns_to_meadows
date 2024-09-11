###############################################################

#From lawns to meadows: Spiders (Arachnida: Araneae) as indicators to 
#measure urban grassland restoration success

#Bach et al. 2024, Urban Ecosystems


###############################################################

# R code to reproduce results

# Author: Alexander Bach
# Last update: 01.09.2024

###############################################################

library(here)
library(data.table)
library(dplyr)
library(vegan)
library(indicspecies)
library(FD)
library(ggplot2)
library(ggpubr)


#### load and prepare necessary data ####
data <- read.csv2(here("spider_data.csv"), sep = "\t", encoding = "UTF-8") 

sites <- read.csv2(here("sites.csv"), sep = "\t", encoding = "UTF-8") 
sites <- sites[ order(sites$LocalityDescription), ]

traits <- read.csv2(here("traits.csv"), sep = "\t", encoding = "UTF-8", row.names = 1, stringsAsFactors = TRUE)

# sum up sampledata from different Collection timepoints
data <- group_by(data, LastIdentificationCache, LocalityDescription_2, CollectionTimeSpan, Value)
data <- summarise(data, NumberOfUnits = sum(NumberOfUnits))

# normalize activity density data according to Saska et al. (2021)
data$NumberOfUnits <- ((data$NumberOfUnits/data$CollectionTimeSpan)/data$Value)

# create crosstable for analyses
data_cr <- setDT(data) %>%
  dcast(LocalityDescription_2~LastIdentificationCache, value.var = "NumberOfUnits", sum)

data_cr <- as.data.frame(data_cr)
row.names(data_cr) <- data_cr$LocalityDescription_2
data_cr$LocalityDescription_2 <- NULL

data_cr <- data_cr[ order(row.names(data_cr)), ]

#### indicator species analysis (ISA) ####

# first exclusive ISA 
IndVal <- multipatt(data_cr, sites$Naturalness, duleg = TRUE, control = how(nperm=999))

# preparation of results
indicators <- as.data.frame(IndVal$sign)
indicators$species = rownames(indicators) 
rownames(indicators)=NULL 
indicators <- indicators[order(indicators$index), ] 
colnames(indicators) <- gsub('^s\\.', "", colnames(indicators))

names(indicators)[(length(colnames(indicators))-3):(length(colnames(indicators)))] <- c("group index", "association_value", "pvalue", "species")

iv <- (indicators$association_value)^2 
indicators <- cbind(indicators, iv)
Naturalness <- indicators$`group index` 
Naturalness <- colnames(indicators)[Naturalness]

indicators <- as.data.frame(cbind(indicators, Naturalness)) 
names(indicators) [ncol(indicators)] <- "Naturalness"

# filter for significance
indicators_sig <- subset(indicators, indicators$iv >= 0.75 & indicators$pvalue<=0.05, check.names = F) 
indicators_sig <- indicators_sig[, -(1:(ncol(indicators_sig)-4))]

# exclude results from dataset
data_cr_ind <- select(data_cr, -one_of(indicators_sig$species))

# second ISA with combined site groups urban vs non urban
IndVal <- multipatt(data_cr_ind, sites$Urban, duleg = TRUE, control = how(nperm=999))

indicators2 <- as.data.frame(IndVal$sign)
indicators2$species = rownames(indicators2) 
rownames(indicators2)=NULL 
indicators2 <- indicators2[order(indicators2$index), ] 
colnames(indicators2) <- gsub('^s\\.', "", colnames(indicators2))

names(indicators2)[(length(colnames(indicators2))-3):(length(colnames(indicators2)))] <- c("group index", "association_value", "pvalue", "species")

iv <- (indicators2$association_value)^2 
indicators2 <- cbind(indicators2, iv)
Urban <- indicators2$`group index` 
Urban <- colnames(indicators2)[Urban]

indicators2 <- as.data.frame(cbind(indicators2, Urban)) 
names(indicators2) [ncol(indicators2)] <- "Urban"

# filter for significance
indicators_sig_urb <- subset(indicators2, indicators2$iv >= 0.75 & indicators2$pvalue<=0.05, check.names = F) 
indicators_sig_urb <- indicators_sig_urb[, -(1:(ncol(indicators_sig_urb)-4))]

# exclude results from dataset
data_cr_ind <- select(data_cr_ind, -one_of(indicators_sig_urb$species))

# third ISA with combined site groups lawn vs non lawn
IndVal <- multipatt(data_cr_ind, sites$Lawn, duleg = TRUE, control = how(nperm=999))

indicators3 <- as.data.frame(IndVal$sign)
indicators3$species = rownames(indicators3) 
rownames(indicators3)=NULL 
indicators3 <- indicators3[order(indicators3$index), ] 
colnames(indicators3) <- gsub('^s\\.', "", colnames(indicators3))

names(indicators3)[(length(colnames(indicators3))-3):(length(colnames(indicators3)))] <- c("group index", "association_value", "pvalue", "species")

iv <- (indicators3$association_value)^2 
indicators3 <- cbind(indicators3, iv)
Lawn <- indicators3$`group index` 
Lawn <- colnames(indicators3)[Lawn]

indicators3 <- as.data.frame(cbind(indicators3, Lawn)) 
names(indicators3) [ncol(indicators3)] <- "Lawn"

# filter for significance
indicators_sig_lawn <- subset(indicators3, indicators3$iv >= 0.75 & indicators3$pvalue<=0.05, check.names = F) 
indicators_sig_lawn <- indicators_sig_lawn[, -(1:(ncol(indicators_sig_lawn)-4))]

# create complete list of indicators
indicators_final <- bind_rows(indicators_sig, indicators_sig_urb, indicators_sig_lawn)

write.csv2(indicators_final, here("indicators_final.csv"))


#### correspondence analysis ####

# data preparation
sites$Naturalness <- as.factor(sites$Naturalness)

# ca
CA_Neu <- cca(data_cr, scale = FALSE)

tiff(filename = "CA_plot.tiff", width = 8, height = 7, units = "in", res = 300)

# plot options
colvec <- c("tan4", "deepskyblue", "darkmagenta", "orange")
cex <- 1.3
par(cex.lab=cex, cex.axis=cex, cex.main=cex)
plot(c(-3.5, 2), c(-4, 2), xlab="CA1 (20.84 %)", ylab="CA2 (16.27%)", type="n") 

ordiellipse(CA_Neu, sites$Naturalness, col=colvec, lwd = 1, lty = 0, draw = "polygon", kind = "ehull")


with(sites, legend("topright", legend = levels(Naturalness), bty = "n", cex = 1.3, pt.cex = 2.5,
                   col = colvec, pch = 19, pt.bg = colvec))

species_score <- scores(CA_Neu)
species_score <- as.data.frame(species_score$species)

subset_spec_score <- subset(species_score, rownames(species_score) %in% indicators_final$species)
subset_spec_score$spec <- row.names(subset_spec_score) %>% make.cepnames()
row.names(subset_spec_score) <- subset_spec_score$spec
subset_spec_score$spec <- NULL

 #add points for indicator species
with(subset_spec_score, points(subset_spec_score, col = "black", pch = 19, cex = 0.8))

# add text labels for indicator species
text(subset_spec_score, labels = rownames(subset_spec_score), col = "black", cex = 1.3, adj = c(-0.3, 0.5))

dev.off()


#### Diversity Measures ####

# calucalte biodiversity parameters
LocalityDescription <- rownames(data_cr)
Shannon <- diversity(data_cr)
Inv <- diversity(data_cr, "inv")
Simpsons <- diversity(data_cr, "simpson")

SpeciesRichness <- specnumber(data_cr)
Abundance <- rowSums(data_cr)


# create dataframe
spec_div <- data.frame(LocalityDescription, SpeciesRichness, Abundance, Simpsons, Shannon)
spec_div <- merge(spec_div, sites, by = "LocalityDescription")
spec_div$Naturalness <- factor(spec_div$Naturalness, levels = c("UI", "UE", "AM", "HM"))
spec_div_melt <- melt(spec_div)

#plot
ggplot(spec_div_melt, aes(x = Naturalness, y = value)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +
  scale_x_discrete(limits = c("UI", "UE", "AM", "HM"))


comp_abu <- list( c("UI", "UE"), c("UE", "HM"), c("UI", "HM"))
comp_div <- list( c("UI", "UE"), c("UE", "HM"), c("UI", "HM"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

abu <- ggboxplot(spec_div, x = "Naturalness", y = "Abundance", fill = "lightgrey", add = "mean_se", ylab = "Daily activity density", xlab = FALSE) + 
  font("x.text", size = 15, color = "black", face = "bold") +
  font("y.text", size = 15, color = "black", face = "bold") +
  font("y.title", size = 15, color = "black", face = "bold") +
  stat_compare_means(comparisons = comp_abu, symnum.args = symnum.args) 

div <- ggboxplot(spec_div, x = "Naturalness", y = "SpeciesRichness", fill = "lightgrey", add = "mean_sd", ylab = "Species richness", xlab = FALSE) + 
  font("x.text", size = 15, color = "black", face = "bold") +
  font("y.text", size = 15, color = "black", face = "bold") +
  font("y.title", size = 15, color = "black", face = "bold") +
  stat_compare_means(comparisons = comp_div, symnum.args = symnum.args) 

ggarrange(abu, div,
          labels = c("A", "B"),
          ncol = 2, nrow = 2)

#### body size distribution analysis ####

data$LocalityDescription <- data$LocalityDescription_2
data$LocalityDescription_2 <- NULL
bs_density <- merge(data, sites, by = "LocalityDescription", all.x = TRUE)
traits$LastIdentificationCache <- row.names(traits)
bs_density <- merge(bs_density, traits, by = "LastIdentificationCache", all.x = TRUE)
traits$LastIdentificationCache <- NULL

# Summarize data by grouping on body.size.female and summarizing NumberOfUnits
bs_summarized <- bs_density %>%
  group_by(Naturalness, body.size.female) %>%
  summarize(NumberOfUnits = sum(NumberOfUnits))

bs_filtered <- bs_summarized %>%
  filter(Naturalness %in% c("UI", "HM", "UE"))

label <- "A: Linyphiidae\nB: Lycosidae (Pardosa sp.)\nC: Alopecosa pulverulenta\nD: Lycosidae (Trochosa sp.)"
colvec_d <- c("deepskyblue", "darkmagenta", "orange")

p <- ggplot(bs_filtered, aes(x = body.size.female, fill = Naturalness, xmin = 1.35)) +
      stat_density(aes(weight = NumberOfUnits), kernel = "gaussian", geom = "density", adjust = 0.5) +
      scale_fill_manual(values=colvec_d) +
      scale_y_continuous(expand = c(0,0)) +
      #geom_vline(xintercept = 1.5, linetype = "dashed") +
      labs(y= "Density", x = "Body Size (female) in mm") +
      theme_pubr(legend= c(0.65,0.8)) +
      font("x.title", size = 15, color = "black", face = "bold") +
      font("y.title", size = 15, color = "black", face = "bold") +
      font("y.text", size = 15, color = "black", face = "bold") +
      font("x.text", size = 15, color = "black", face = "bold") +
      geom_rect(aes(xmin = 1.3, xmax = 3.3, ymin = 0.3, ymax = 0.52),
            fill = "transparent", color = "red", linewidth = 1.3) +
      geom_rect(aes(xmin = 5, xmax = 7.5, ymin = 0.2, ymax = 0.58),
           fill = "transparent", color = "red", linewidth = 1.3) +
      geom_rect(aes(xmin = 8.5, xmax = 9.8, ymin = 0.08, ymax = 0.15),
            fill = "transparent", color = "red", linewidth = 1.3) +
      geom_rect(aes(xmin = 10.3, xmax = 11.8, ymin = 0.02, ymax = 0.13),
            fill = "transparent", color = "red", linewidth = 1.3)

p <- ggpar(p,font.legend = c(15, "bold", "black"), legend.title = " ", xticks.by = 2, xlim = c(1.3,15)) +
        annotate("text", x = 1.6, y = 0.5, label = "A", size = 8) +
        annotate("text", x = 5.3, y = 0.56, label = "B", size = 8) +
        annotate("text", x = 8.75, y = 0.135, label = "C", size = 8) +
        annotate("text", x = 10.56, y = 0.115, label = "D", size = 8) +
        annotate("text", x = 10, y = 0.32, label = label, size = 5, hjust = 0, fontface = 1) +
        annotate("text", x = 10, y = 0.38, label = paste("Species or Groups causing density peaks"), size = 5, 
           hjust = 0, fontface = 2)

ggsave("my_plot_cor.tiff", width = 2306, height = 1446, units = "px", plot = p, device = "tiff")
pdf("my_plot_cor.pdf", width = 12.01, height = 7.531)
plot(p)
dev.off()

#### traits analyses ####
traits$dispersal <- as.numeric(traits$dispersal)
CWM_traits <- functcomp(traits, as.matrix(data_cr), CWM.type = "all")
CWM_traits$LocalityDescription <- rownames(CWM_traits)
CWM_traits_melt <- reshape2::melt(CWM_traits, id = "LocalityDescription")
CWM_traits_melt <- base::merge(CWM_traits_melt, sites, by = "LocalityDescription") #add Site Properties for Grouping

CWM_traits <- cbind(CWM_traits, sites$Naturalness)
CWM_traits$Naturalness <- CWM_traits$`sites$Naturalness`
CWM_traits$LocalityDescription <- NULL
CWM_traits$`sites$Naturalness` <- NULL
CWM_traits$Naturalness <- factor(CWM_traits$Naturalness, levels = c("UI", "UE", "AM", "HM"))

# loop through each column in CWM_traits (excluding the last column, which is the grouping variable) with anova and Tukey HSD
for (i in 1:(ncol(CWM_traits)-1)) {
  
  if(i == 1){
    results <- data.frame()
    dif_hsd <- array(1:6, dim = 6) 
  }
  
  # get the name of the current trait
  trait <- colnames(CWM_traits)[i]
  
  # perform a ANOVA to compare the trait values in the UI group to the rest of the groups
  sig_test <- aov(CWM_traits[,i] ~ Naturalness, data = CWM_traits)
  sum_test = unlist(summary(sig_test))
  p_val <- sum_test["Pr(>F)1"]
  f_val <- sum_test["F value1"]
  
  if (p_val < 0.05) {
    dif <- "Different"
  } else {
    dif <- "Not Different"
  }
  
  # perform pairwise comparisons using TukeyHSD
  pairwise_comparisons <- TukeyHSD(sig_test)
  pairwise_comparisons <- as.data.frame(pairwise_comparisons$Naturalness)
  pairwise_comparisons$groups <- row.names(pairwise_comparisons)
  
  # extract the relevant information from TukeyHSD 
  group <- pairwise_comparisons$groups
  p_val_hsd <- pairwise_comparisons$`p adj`
  val_size <- pairwise_comparisons$diff
  
  for(j in 1:length(p_val_hsd)) {
  if (p_val_hsd[j] < 0.05) {
    dif_hsd[j] <- "Different"
  } else {
    dif_hsd[j] <- "Not Different"
  }}
  
  results <- rbind(results, data.frame(p_val = rep(p_val, length(pairwise_comparisons$groups)),
                                       f_val = rep(f_val, length(pairwise_comparisons$groups)),
                                       dif = rep(dif,length(pairwise_comparisons$groups)),
                                       trait = rep(trait, length(pairwise_comparisons$groups)),
                                       p_hsd = p_val_hsd,
                                       dif_hsd = dif_hsd,
                                       group = pairwise_comparisons$groups,
                                       effect = val_size))
  

  if(i == (ncol(CWM_traits)-1)){
    sig_results <- subset(results, dif_hsd == "Different" & grepl("UE", group) & !grepl("UI", group))
    sig_results$Winner <- ifelse(sig_results$effect > 0, substr(sig_results$group, 1, 2), substr(sig_results$group, 4, 5))
  }
  
}