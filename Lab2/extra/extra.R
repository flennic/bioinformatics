# Lab2 - Assignment 1 - question 1: R script of the extra documentation provided
# NOTE: the final workspace is saved in the file "data_lizard.RData"

# Based on http://www.jcsantosresearch.org/Class_2014_Spring_Comparative/pdf/week_2/Jan_13_15_2015_GenBank_part_2.pdf

seq_1_DNAbin <- read.GenBank("JF806202")
seq_1_DNAbin
class(seq_1_DNAbin)
attr(seq_1_DNAbin, "species") # To get the specie name of the sequence
seq_1_DNAbin$JF806202
str(seq_1_DNAbin) # We get the structure of the object

# Save as character object:
seq_1_character <- read.GenBank("JF806202", as.character = TRUE)
seq_1_character # This is not a very nice format

lizards_accession_numbers <- c("JF806202", "HM161150", "FJ356743", "JF806205",
                               "JQ073190", "GU457971", "FJ356741", "JF806207",
                               "JF806210", "AY662592", "AY662591", "FJ356748",
                               "JN112660", "AY662594", "JN112661", "HQ876437",
                               "HQ876434", "AY662590", "FJ356740", "JF806214",
                               "JQ073188", "FJ356749", "JQ073189", "JF806216",
                               "AY662598", "JN112653", "JF806204", "FJ356747",
                               "FJ356744", "HQ876440", "JN112651", "JF806215",
                               "JF806209")   # Create a vector a GenBank accession numbers

# Read sequences and place them in a DNAbin object
lizards_sequences <- read.GenBank(lizards_accession_numbers) 
lizards_sequences # A brief summary of what is in the object, including base composition
str(lizards_sequences) # A list of the DNAbin elements with length of the sequences
# Notice the one of the attributes is the species names

attributes(lizards_sequences) # See the list of attributes and contents
names(lizards_sequences) # The accession numbers
attr(lizards_sequences, "species") # We get the species list. 
# Notice this attr(.) is slightly different function

lizards_sequences_GenBank_IDs <- paste(attr(lizards_sequences, "species"), 
                                       names(lizards_sequences), sep ="_RAG1_")
## build a character vector with the species, GenBank accession numbers, and gene
## name "_RAG1_” this is its common abbreviation: recombination activating protein 1
## notice the use of the paste function: textA, textB, textC
## results in: textAtextCtextB
lizards_sequences_GenBank_IDs  # A more informative vector of names for our sequences

?write.dna # This function writes in a file a list of DNA sequences in sequential,
# interleaved, or FASTA format.
### we are going to write in fasta format
write.dna(lizards_sequences, file ="lizard_fasta_1.fasta", format = "fasta", append =
            FALSE, nbcol = 6, colsep = " ", colw = 10)
########### Some relevant arguments for write.dna()
#x: a list or a matrix of DNA sequences.
#file: a file name specified to contain our sequences
#format: Three choices are possible: "interleaved", "sequential", or "fasta", or any
#unambiguous abbreviation of these.
#append: a logical, if TRUE the data are appended to the file without erasing the data
#possibly existing in the file, otherwise the file is overwritten (FALSE the default).
#nbcol: a numeric specifying the number of columns per row (6 by default)
#colsep: a character used to separate the columns (a single space by default).
#colw: a numeric specifying the number of nucleotides per column (10 by default).
###########

# Read our fasta file using the seqinr package
lizard_seq_seqinr_format <- read.fasta(file = "lizard_fasta_1.fasta", seqtype = "DNA",
                                       as.string = TRUE, forceDNAtolower = FALSE)
lizard_seq_seqinr_format   # This shows different form to display the same sequence information 

write.fasta(sequences = lizard_seq_seqinr_format, names = lizards_sequences_GenBank_IDs,
            nbchar = 10, file.out = "lizard_seq_seqinr_format.fasta")
# Suggestion: Do not rearrange, delete or add sequenced to the fasta file, as the
# function will assign the names in the order provided in the file and the name vector

# We can use a package that use an API (application programming interface) 
# to interact with the NCBI website. More info in: http://en.wikipedia.org/wiki/Application_programming_interface
# RUN ONLY ONCE: install.packages ("rentrez")
library(rentrez)

lizard <- "Basiliscus basiliscus[Organism]" # We want a character vector

# Nucleotide database (nuccore) and retmax determines the max number
lizard_search <- entrez_search(db="nuccore", term=lizard, retmax=40)
lizard_search
lizard_search$ids # Gives you the NCBI ids
# Gets your sequences as a character vector
# NOTE: you may need to run the following line in advance:
# httr::set_config(httr::config(http_version = 0))
lizard_seqs <- entrez_fetch(db="nuccore", id=lizard_search$ids, rettype="fasta")
lizard_seqs

Bbasiliscus_RAG1 <- "Basiliscus basiliscus[Organism] AND RAG1[Gene]"
Bbasiliscus_RAG1_search <- entrez_search(db="nuccore", term=Bbasiliscus_RAG1, retmax=10)
# Nucleotide database (nuccore) and retmax determines no more than 10 access numbers to return
Bbasiliscus_RAG1_search$ids #gives you the NCBI ids
Bbasiliscus_RAG1_seqs <- entrez_fetch(db="nuccore", id=Bbasiliscus_RAG1_search$ids,
                                      rettype="fasta")
Bbasiliscus_RAG1_seqs # Notice \n (new line) delimiter. Other common delimiters are \r
# (carriage return) and \t (tab).
write(Bbasiliscus_RAG1_seqs, "Bbasiliscus_RAG1.fasta", sep="\n") # Gets sequence to a file

Bbasiliscus_RAG1_seqinr_format <- read.fasta(file = "Bbasiliscus_RAG1.fasta", seqtype =
                                               "DNA", as.string = TRUE, forceDNAtolower = FALSE)
Bbasiliscus_RAG1_seqinr_format # you can also check the .fasta file in the working folder

# We can use the ‘rentrez’ package to get lots of sequences using taxonomic 
# classifications for specific markers
Liolaemus_CYTB <- "Liolaemus[Organism] AND CYTB[Gene]"
# This is a well-studied gene from this genus of South American lizards
Liolaemus_CYTB_search <- entrez_search(db="nuccore", term=Liolaemus_CYTB, retmax=100)
Liolaemus_CYTB_search # There are 2539 sequences that match this query

Liolaemus_CYTB_search_2 <- entrez_search(db="nuccore", term=Liolaemus_CYTB, retmax=2539)
Liolaemus_CYTB_search_2$ids # Gives you the NCBI ids
Liolaemus_CYTB_seqs <- entrez_fetch(db="nuccore", id=Liolaemus_CYTB_search_2$ids ,
                                    rettype="fasta")
# We get an error “client error: (414) Request-URI Too Long”. We are asking too many sequences

# Lets adjust the search and fetch by smaller chunks so we can get the first 1500 sequences
# NOTE: in the original PDF the split was every 500 sequences. 
# It was too long and it was necessary to split it.
Liolaemus_CYTB_seqs_part_1 <- entrez_fetch(db="nuccore", id=Liolaemus_CYTB_search_2$ids
                                           [1:300] , rettype="fasta")
Liolaemus_CYTB_seqs_part_2 <- entrez_fetch(db="nuccore", id=Liolaemus_CYTB_search_2$ids
                                           [301:600] , rettype="fasta")
Liolaemus_CYTB_seqs_part_3 <- entrez_fetch(db="nuccore", id=Liolaemus_CYTB_search_2$ids
                                           [601:900] , rettype="fasta")
Liolaemus_CYTB_seqs_part_4 <- entrez_fetch(db="nuccore", id=Liolaemus_CYTB_search_2$ids
                                           [901:1200] , rettype="fasta")
Liolaemus_CYTB_seqs_part_5 <- entrez_fetch(db="nuccore", id=Liolaemus_CYTB_search_2$ids
                                           [1201:1500] , rettype="fasta")

write(Liolaemus_CYTB_seqs_part_1, "Liolaemus_CYTB_seqs.fasta", sep="\n")
write(Liolaemus_CYTB_seqs_part_2, "Liolaemus_CYTB_seqs.fasta", sep="\n", append = TRUE)
# It gets the sequences to the same file by changing the logical argument of append from
# the default FALSE to TRUE (i.e., can abbreviate TRUE with T or other unambiguous abbreviation)
write(Liolaemus_CYTB_seqs_part_3, "Liolaemus_CYTB_seqs.fasta", sep="\n", append = TRUE)
write(Liolaemus_CYTB_seqs_part_4, "Liolaemus_CYTB_seqs.fasta", sep="\n", append = TRUE)
write(Liolaemus_CYTB_seqs_part_5, "Liolaemus_CYTB_seqs.fasta", sep="\n", append = TRUE)
# You will get a 1.3 Mb file with all 1500 sequences

# We can read our fasta file using the seqinr package and rename the sequences
Liolaemus_CYTB_seqs_seqinr_format <- read.fasta(file = "Liolaemus_CYTB_seqs.fasta",
                                                seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
Liolaemus_CYTB_seqs_seqinr_format
Liolaemnus_CYTB_names <- attr(Liolaemus_CYTB_seqs_seqinr_format, "name")
Liolaemnus_CYTB_names
Liolaemnus_CYTB_names <- gsub("\\..*","", Liolaemnus_CYTB_names)
# Eliminate characters after "." using ?gsub (Pattern Matching and Replacement)
Liolaemnus_CYTB_names <- gsub("^.*\\|", "", Liolaemnus_CYTB_names)
# Eliminate characters before "|" using ?gsub (Pattern Matching and Replacement)
Liolaemnus_CYTB_names

# We can read our fasta file using ape package to get accession numbers and species names
# NOTE: it takes time
Liolaemus_CYTB_seqs_ape_format <- read.GenBank(Liolaemnus_CYTB_names)
attr(Liolaemus_CYTB_seqs_ape_format, "species")
# To get the species names of the sequence
names(Liolaemus_CYTB_seqs_ape_format)
Liolaemus_CYTB_seqs_GenBank_IDs <- paste(attr(Liolaemus_CYTB_seqs_ape_format,
                                              "species"), names(Liolaemus_CYTB_seqs_ape_format), sep="_CYTB_")
## build a vector object with the species, GenBank accession numbers, and type of gene
Liolaemus_CYTB_seqs_GenBank_IDs # Vector of names to add to sequences
# Read our fasta file 'Liolaemus_CYTB_seqs.fasta' using seqinr package
Liolaemus_CYTB_seqs_seqinr_format <- read.fasta(file = "Liolaemus_CYTB_seqs.fasta",
                                                seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
# Rewrite our fasta file using the name vector that we created previously
write.fasta(sequences = Liolaemus_CYTB_seqs_seqinr_format, names =
              Liolaemus_CYTB_seqs_GenBank_IDs, nbchar = 10, file.out =
              "Liolaemus_CYTB_seqs_seqinr_format.fasta")




######################################################################################################

# Lab2 - Assignment 1 - question 2: R script of the extra documentation (first link) suggested

# Create a tree
tt <- "(Lemur,(Tarsius,(((Callithrix,Saimiri),(Pithecia,Lagothrix)),((Macaca,Colobus),
       (Hylobates,(Pongo,(Gorilla,(Pan,Homo))))))));"
pr.tree <- read.tree(text = tt)
# You can see individual components of this list by typing:
pr.tree$edge    # or
pr.tree$tip.label

write.tree(pr.tree)
write.nexus(pr.tree, file = "primate.tre")
