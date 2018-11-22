```{r fig.width = 8, fig.height=8}

library(phangorn)

# ------------------------------------------------------------------------------
# Question 1.2
# ------------------------------------------------------------------------------

# Simulate phylogenetic tree with 33 tips in phylo format (ape) ----------------
set.seed(1)
tree = ape::rtree(n = 33)

# Plot resulting tree ----------------------------------------------------------

plot(tree, edge.width = 1, main = "Plot of simulated phylogenetic tree")
# phytools::plotTree(tree) # Alternative

# Simulate sequences on this tree using phangorn::simSeq() ---------------------
Q = matrix(c(.1, .8, .05, .05, 
             .35, .1, .1, .45,
             .3, .2, .2, .3,
             .6, .1, .25, .05), nrow = 4, byrow = TRUE)
rownames(Q) = c("a", "c", "g", "t")
colnames(Q) = c("a", "c", "g", "t")

tree_sequences_sim = phangorn::simSeq(tree, l = 2000, Q = Q, bf = Original)

# Explanation of parameters: 
# l = 2000 because average sequence length in given data is ca. 2000 
# bf = Original because this is the vector with the original base proportions
# Q = just chosen the matrix from Special Exercise 1 (Question 3)

# Convert to DNAbin
tree_sequences_sim = as.DNAbin(tree_sequences_sim)

# Save simulated sequences as fasta file ---------------------------------------

# Write simulated lizard sequences as fasta file
ape::write.dna(tree_sequences_sim, file ="lizard_seqs_tree_sim.fasta", format = "fasta", 
               append = FALSE, nbcol = 6, colsep = " ", colw = 10)

# Report base composition ------------------------------------------------------

# Check base distribution. Note: very similar but different (as expected).
Original = ape::base.freq(lizards_sequences)
Simulated_Tree = ape::base.freq(tree_sequences_sim)
knitr::kable(data.frame(Original, Simulated_Tree), 
             caption = "Base composition")

```