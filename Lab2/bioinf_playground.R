library(msa)
### Set working directory!
setwd("./bioinformatics")
### Install 'Biostrings' and 'msa' packages

# Import the fasta files
original_as_AAStringSet =
  readAAStringSet("data/lizard_seqs.fasta", format = "fasta")
sim_sample_as_AAStringSet =
  readAAStringSet("data/simulated_lizards.fasta", format = "fasta")
sim_tree_as_AAStringSet =
  readAAStringSet("data/lizard_seqs_tree_sim.fasta", format = "fasta")

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
                compress = FALSE, compression_level = NA, format = "fasta", ...)

writeXStringSet(as(sim_tree_aligned, "XStringSet"),
                "data/simulated_lizards_tree_aligned.fasta", append = FALSE,
                compress = FALSE, compression_level = NA, format = "fasta", ...)
