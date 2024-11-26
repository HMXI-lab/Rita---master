# Why this code?
## To process and analyse Amplicon Sequence Variants (ASVs) from microbiome data, particularly focusing on 16S rRNA gene sequences 

# Load necessary libraries
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install(c("Biostrings", "msa", "ape", "seqinr"))


# # Install and load necessary packages
# if (!requireNamespace("rentrez", quietly = TRUE)) {
#   install.packages("rentrez")
# }


library(Biostrings) #For handling biological sequences
library(msa) # For running Clustal Omega or MAFFT (multiple sequence alignment)
library(ape) # For distance calculation = phylogenetic analysis
library(phangorn)  # For efficient distance calculations, through phylogenetic tree construction and, after, distance calculations
library(igraph)    # For graph-based clustering
library(seqinr)    # For sequence retrieval and analysis
library(tidyverse)  # For data manipulation and visualization
library(ggplot2)    # For data manipulation and visualization



setwd("~/ITQB/Omics")

#To remove ASVs (= Amplicon Sequence Variant) likely due to contamination, filter the ASVs that have higher than 1000 reads as average across samples

ASVs_Rita <- read_tsv("Test_run_ASV_table_240_240.tsv")
View(ASVs_Rita) 

# Assuming your data frame is called df
df_filtered <- ASVs_Rita %>%
  # Calculate the sum of counts for each ASV column (excluding the first column)
  summarise(across(-Sample, mean)) %>%
  # Keep only the columns where the sum is greater than 900
  select(where(~ . > 1000)) 

# Filter the original data frame to keep the sample_id and filtered ASV columns
df_filtered <- ASVs_Rita %>% 
  select(Sample, names(df_filtered))

#Only 12 ASVs are kept after the filtering, then, we select only 12 ASVs from the taxonomy table, and look for ASVs that could
#have similarity > 98% (do clustalW, calculate similarity matrix, and report those ASVs with sequences higher than 98% identity)


# Step 1: Load the taxonomy table
taxonomy_table <- read_tsv("Test_run_ASV_taxonomy_240_240.tsv", col_types = cols(.default = "c"))

# Keep only ASVs 1 to 15
taxonomy_table <- taxonomy_table[1:12, ]

# Step 2: Create a DNAStringSet object with named sequences
sequences <- DNAStringSet(taxonomy_table$Seq)
names(sequences) <- taxonomy_table$ID  # Assign ASV IDs as names

# Step 3: Perform ClustalW alignment using the msa package
aligned <- msaClustalW(sequences)

# Convert aligned sequences to a character vector and keep names
aligned_sequences <- as.character(aligned)
aligned_names <- names(aligned_sequences)

# Step 4: Compute pairwise sequence similarity with correct ASV IDs
compute_pairwise_similarity <- function(aligned_sequences, aligned_names) {
  num_seqs <- length(aligned_sequences)
  similarity_matrix <- matrix(0, nrow = num_seqs, ncol = num_seqs)
  rownames(similarity_matrix) <- aligned_names
  colnames(similarity_matrix) <- aligned_names
  
  for (i in 1:num_seqs) {
    for (j in i:num_seqs) {
      seq1 <- unlist(strsplit(aligned_sequences[i], ""))
      seq2 <- unlist(strsplit(aligned_sequences[j], ""))
      matches <- sum(seq1 == seq2)
      identity <- (matches / length(seq1)) * 100
      similarity_matrix[i, j] <- identity
      similarity_matrix[j, i] <- identity
    }
  }
  return(similarity_matrix)
}

similarity_matrix <- compute_pairwise_similarity(aligned_sequences, aligned_names)
similarity_matrix


# Step 5: Identify clusters of similar ASVs and prepare the report
threshold <- 98.0  # Threshold for collapsing
collapse_report <- data.frame(ASV_ID_1 = character(), ASV_ID_2 = character(), Percent_Identity = numeric())

for (i in 1:(nrow(similarity_matrix) - 1)) {
  for (j in (i + 1):ncol(similarity_matrix)) {
    if (similarity_matrix[i, j] >= threshold) {
      collapse_report <- rbind(collapse_report, data.frame(
        ASV_ID_1 = rownames(similarity_matrix)[i],
        ASV_ID_2 = colnames(similarity_matrix)[j],
        Percent_Identity = similarity_matrix[i, j]
      ))
    }
  }
}

# Save the report to a CSV file
write.csv(collapse_report, "collapse_report.csv", row.names = FALSE)

print("Collapse report generated successfully as 'collapse_report.csv'.")


# We obtained two ASVs with similarity between them higher than 98%: ASV0006 and ASV0010
# Let's do nucleotide blast to confirm whether they align with the same or different species

# Example ASV IDs with high similarity (replace these with your actual ASV IDs)
high_similarity_asv1 <- "ASV0006"  # Replace with the actual ASV ID
high_similarity_asv2 <- "ASV0010"  # Replace with the actual ASV ID

# Extract sequences corresponding to the ASVs with high similarity
sequence1 <- taxonomy_table$Seq[taxonomy_table$ID == high_similarity_asv1]
sequence2 <- taxonomy_table$Seq[taxonomy_table$ID == high_similarity_asv2]

#copy each sequence and do a blastn in NCBI
sequence1
sequence2


#In this case ASV006 or sequence1 is Bacteroides fragilis;
#ASV0010 aligns with Bacteroides thetaiotaomicron


#As no sequences are to be collapsed, do the merged table


ASVs_Rita_v1 <-gather(data = df_filtered, key = "ID", value = "counts", -Sample)
ASVs_Rita_v1


ASVs_Tax_Rita <- left_join(ASVs_Rita_v1, taxonomy_table)
ASVs_Tax_Rita

#Add the species names mannually

manual_updates <- data.frame(
      ID = c("ASV0001",
             "ASV0002",
             "ASV0003",
             "ASV0004",
             "ASV0005",
             "ASV0006",
             "ASV0007",
             "ASV0008",
             "ASV0009",
             "ASV0010",
             "ASV0011",
             "ASV0012"),  # Replace with actual ASV IDs
  Species = c("Erysipelatoclostridium ramosum",
              "Escherichia coli",
              "Streptococcus salivarius",
              "Enterococcus faecalis",
              "Lactobacillus fermentum",
              "Bacteroides fragilis",
              "Dorea formicigenerans",
              "Coprococcus comes",
              "Clostridium perfringens",
              "Bacteroides thetaiotaomicron",
              "Enterocloster bolteae",
              "Meditirraneibacter gnavus"),  # Replace with actual species names
  stringsAsFactors = FALSE
)


ASVs_Tax_Rita  <- ASVs_Tax_Rita  %>%
  left_join(manual_updates, by = "ID", suffix = c("", ".new")) %>%
  mutate(Species = coalesce(Species.new, Species)) %>%
  select(-Species.new)  # Remove the temporary column



ASVs_Tax_Rita <- ASVs_Tax_Rita %>%
  separate(Sample, into = c("Info1", "BiolRepl", "Community", "Concentration", "Info2", "Info3"), sep = "_", fill = "right", remove = FALSE)



write_tsv(ASVs_Tax_Rita, "filtered_ASVs_Tax_Rita.tsv")


# Define the internal standard ASV ID
internal_standard_asv <- "ASV0005"


# Step 1: Calculate total counts per biological replicate, community and concentration, excluding the internal standard
ASVs_Tax_Rita_counts <- ASVs_Tax_Rita %>%
  group_by(BiolRepl, Community, Concentration) %>%
  mutate(total_counts = sum(counts),
         total_counts_excl_internal = sum(counts[ID != internal_standard_asv])) %>%
  filter(total_counts > 999)

# Step 2: Calculate internal standard counts separately
internal_standard_counts <- ASVs_Tax_Rita_counts %>%
  filter(ID == internal_standard_asv) %>%
  select(BiolRepl, Community, Concentration, internal_standard_count = counts)

# Step 3: Merge internal standard counts with the main data and calculate the relative counts considering and not considering the internal standard
ASVs_Tax_Rita_rel_counts <- ASVs_Tax_Rita_counts %>%
  left_join(internal_standard_counts, by = c("BiolRepl", "Community", "Concentration")) %>%
  group_by(BiolRepl, Community, Concentration) %>%
  mutate(rel_counts_1 = counts / total_counts_excl_internal,
         rel_counts_2 = counts / total_counts,
         rel_counts_normalized = rel_counts_1 / internal_standard_count) %>%
  select(-internal_standard_count)  # Remove the temporary column


#the results with the internal standard are very weird, so, I will use the AUC instead to obtain the absolute abundances


#Open the AUCs table (I modified the table a bit, so it can be merged easily with the taxonomy and counts table)
library(readr)
All_AUC <- read_delim("All_AUC.csv", delim = ";", 
                      escape_double = FALSE, trim_ws = TRUE)
View(All_AUC)

#Calculate AUCs average per technical replicate

All_AUC_1 <- All_AUC %>% group_by(Sample, BiolRepl, Community, drug_conc) %>%
  summarise(meanFinalOD = mean(finalOD), meanAUC = mean (AUC), meanMaxOD = mean(maxOD)) %>% ungroup() %>%
  select (- BiolRepl, -Community, -drug_conc) 

#Prepare the taxonomy-counts table and the AUC tables for merging

ASVs_Tax_Rita_rel_counts <- ASVs_Tax_Rita_rel_counts %>% mutate (Sample1 = Sample) %>%
  separate(Sample, into = c("Sample", "Info4"), sep = "_S", fill = "right", remove = FALSE)


#Merge tables
ASVs_Tax_Rita_rel_counts_ODs  <- ASVs_Tax_Rita_rel_counts   %>%
  left_join(All_AUC_1, by = "Sample")

#Multiply the rel_counts by the AUCs
ASVs_Tax_Rita_rel_counts_ODs_norm <- ASVs_Tax_Rita_rel_counts_ODs %>% 
  mutate(rel_counts_normalized_AUCs = rel_counts_1 * meanAUC)


write_tsv(ASVs_Tax_Rita_rel_counts_ODs_norm, "ASVs_Tax_Rita_rel_counts_ODs_norm.tsv")


library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

ASVs_Tax_Rita_rel_counts_ODs_norm$Concentration <- factor (ASVs_Tax_Rita_rel_counts_ODs_norm$Concentration, levels =c("Zn", "0", "24", "30", "40"))

#Plot only the relative counts
# Create scatter plot with faceting
ASVs_Tax_Rita_rel_counts_ODs_norm %>%
  ggplot(aes(x = Species, y = rel_counts_1, color = BiolRepl)) +
  geom_point() +
  facet_grid(Community ~ Concentration) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  

#Plot only the relative counts normalized by the internal standard (does not look great)
# Create scatter plot with faceting
ASVs_Tax_Rita_rel_counts_ODs_norm %>%
  ggplot(aes(x = Species, y = rel_counts_normalized, color = BiolRepl)) +
  geom_point() +
  facet_grid(Community ~ Concentration) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  

# Create scatter plot with faceting (normalized_AUCs)
ASVs_Tax_Rita_rel_counts_ODs_norm %>%
  ggplot(aes(x = Species, y = rel_counts_normalized_AUCs, color = BiolRepl)) +
  geom_point() +
  facet_grid(Community ~ Concentration) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  


#Calculate the mean of the relative counts
ASVs_Tax_Rita_rel_counts_ODs_norm_mean <- ASVs_Tax_Rita_rel_counts_ODs_norm %>% group_by(Community, Concentration, Kingdom, Phylum, Class, Order,Family,Genus,Species) %>%
  summarise(mean_rel_counts_AUC = mean(rel_counts_normalized_AUCs), sd_rel_counts_AUCs = sd (rel_counts_normalized_AUCs) , mean_rel_counts = mean(rel_counts_1, na.rm=TRUE), sd_rel_counts = sd(rel_counts_2), mean_counts = mean(counts), 
            .groups = 'drop')


ASVs_Tax_Rita_rel_counts_ODs_norm_mean$Concentration <- factor (ASVs_Tax_Rita_rel_counts_ODs_norm_mean$Concentration, levels =c("Zn", "0", "24", "30", "40"))

write_tsv(ASVs_Tax_Rita_rel_counts_ODs_norm_mean , "ASVs_Tax_Rita_rel_counts_ODs_norm_mean.tsv")


#Plot the relative abundances, i.e., not normalized with the AUC
# Create a bar plot with geom_col()
ASVs_Tax_Rita_rel_counts_ODs_norm_mean %>% filter (Species != "Lactobacillus fermentum") %>%
  ggplot(aes(x = Species, y = mean_rel_counts, fill = Species)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_rel_counts - sd_rel_counts, 
                    ymax = mean_rel_counts + sd_rel_counts), 
                width = 0.2) +  
  facet_grid(Community ~ Concentration) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = mycolors) +
  labs(x = "Species", y = "Relative mean Counts", title = "Relative mean Counts by Species")

#Plot the relative abundances, i.e., not normalized with the AUC
# Create a stacked bar plot
ASVs_Tax_Rita_rel_counts_ODs_norm_mean  %>%
  filter(Species != "Lactobacillus fermentum") %>%
  ggplot(aes(x = Community, y = mean_rel_counts, fill = Species)) +
  geom_col(position = "stack") +  # Stacked bars
  facet_grid(. ~ Concentration) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = mycolors) +
  labs(x = "Community", y = "Relative mean Counts", title = "Relative mean Counts by Species and Community")

#Plot the absolute abundances, i.e., normalized with the AUC
# Create a bar plot with geom_col()
ASVs_Tax_Rita_rel_counts_ODs_norm_mean %>% filter (Species != "Lactobacillus fermentum") %>%
  ggplot(aes(x = Species, y = mean_rel_counts_AUC, fill = Species)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_rel_counts_AUC - sd_rel_counts_AUCs, 
                    ymax = mean_rel_counts_AUC + sd_rel_counts_AUCs), 
                width = 0.2) +  
  facet_grid(Community ~ Concentration) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = mycolors) +
  labs(x = "Species", y = "Absolute mean Counts", title = "Absolute mean Counts by Species")

#Plot the absolute abundances, i.e., normalized with the AUC
# Create a stacked bar plot
ASVs_Tax_Rita_rel_counts_ODs_norm_mean  %>%
  filter(Species != "Lactobacillus fermentum") %>%
  ggplot(aes(x = Community, y = mean_rel_counts_AUC, fill = Species)) +
  geom_col(position = "stack") +  # Stacked bars
  facet_grid(. ~ Concentration) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = mycolors) +
  labs(x = "Community", y = "Absolute mean Counts", title = "Absolute mean Counts by Species and Community")




