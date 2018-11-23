library(msa)
### Set working directory!
setwd("./bioinformatics")
### Install 'Biostrings' and 'msa' packages
# These uncommented lines can be used to install the BiocManager which enables you to install Biostringss
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#library("BiocManager")
#BiocManager::install("msa")

# Import the fasta files
original_as_AAStringSet =
  readAAStringSet("data/lizard_seqs.fasta", format = "fasta")
sim_sample_as_AAStringSet =
  readAAStringSet("data/simulated_lizards.fasta", format = "fasta")
sim_tree_as_AAStringSet =
  readAAStringSet("data/simulated_lizards_tree.fasta", format = "fasta")

# Calculate the alignment using ClustalW (default, also possible are 'ClustalOmega' and 'Muscle')
original_aligned = msa(original_as_AAStringSet, method = "ClustalW")
sim_sample_aligned = msa(sim_sample_as_AAStringSet, method = "ClustalW")
sim_tree_aligned = msa(sim_tree_as_AAStringSet, method = "ClustalW")

# Save them to files
writeXStringSet(as(original_aligned, "XStringSet"),
                "data/lizard_seqs_aligned.fasta", append = FALSE,
                compress = FALSE, compression_level = NA, format = "fasta")

writeXStringSet(as(sim_sample_aligned, "XStringSet"),
                "data/simulated_lizards_aligned.fasta", append = FALSE,
                compress = FALSE, compression_level = NA, format = "fasta")

writeXStringSet(as(sim_tree_aligned, "XStringSet"),
                "data/simulated_lizards_tree_aligned.fasta", append = FALSE,
                compress = FALSE, compression_level = NA, format = "fasta")

# This actually gives a nice output! Seems like there a a little bit to many
# spaces in between
library(DECIPHER)
original <- readDNAStringSet("data/lizard_seqs.fasta", format = "fasta")
sim_seq <- readDNAStringSet("data/simulated_lizards.fasta", format = "fasta")
tree_seq <- readDNAStringSet("data/simulated_lizards_tree.fasta", format = "fasta")
#original <- OrientNucleotides(seqs)
aligned_original <- AlignSeqs(original)
aligned_sim <- AlignSeqs(sim_seq)
set.seed(1)
tree = ape::rtree(n = 33)
aligned_tree <- AlignSeqs(tree_seq, guideTree = phylogram::as.dendrogram(tree))
BrowseSeqs(aligned_original, highlight=0)
BrowseSeqs(aligned_sim, highlight=0)
BrowseSeqs(aligned_tree, highlight=0)

# Calculate distance matrix (no option to pick the type of measure, it is only dissimilarity)
dist_original <- DistanceMatrix(original)
dist_sim <- DistanceMatrix(sim_seq)
dist_tree <- DistanceMatrix(tree_seq)

Heatmap(dist_original, column_title = "Original Sequences",
        show_row_names = FALSE)
Heatmap(dist_sim, column_title = "Simulated Sequences",
        show_row_names = FALSE)
Heatmap(dist_tree, column_title = "Simulated Sequences (Tree)",
        show_row_names = FALSE)
