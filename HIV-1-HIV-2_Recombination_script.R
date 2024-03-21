setwd("D:/Dropbox/PhD 2023 UG/In silico/Recombination project/nearful_HIV1_2/phylogenetic tree of potential recombinants/sieme_script")

########## install necessary packges 
install.packages("remotes")
remotes::install_github("GuangchuangYu/treeio")
install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("DECIPHER")
install.packages("ape")
install.packages("adegenet")
install.packages("Biostrings")
install.packages("ggtree")
install.packages("treeio")
install.packages("pheatmap")

#######Load library############
library(seqinr)
library(adegenet)
library(ape)
library(ggplot2)
library(ggtree)
library(DECIPHER)
library(viridis)
library(pheatmap)
library(treeio)
library(Biostrings)
library(pheatmap)

# load the sequences from a text file

seq <- readDNAStringSet("potential_rec.txt", format = "fasta") #specify file format as fasta
seqs <- DNAStringSet(seq) ###Create a DNAStringSet object

# look at some of the sequences (optional)
seqs

print(seqs, show = "complete") #look at the summary characteristics of your sequences

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform the alignment using DECIPHER
aligned <- AlignSeqs(seqs)
aligned #view alignment

# view the alignment in a browser (optional) using DECIPHER 
BrowseSeqs(aligned, highlight=0)

# write the alignment to a new FASTA file using Biostrings
writeXStringSet(aligned,
                file="potential_rec_aligned.fasta")

# read in the aligned data using seqinr 
dna <- read.alignment("potential_rec_aligned.fasta", format = "fasta")
dna

# create a distance matrix for the alignment 
distMat <- dist.alignment(dna, matrix = "similarity")
distMat

temp <- as.data.frame(as.matrix(distMat))
temp

# Customized heatmap using pheatmap
# Assuming you have the 'pheatmap' package installed; otherwise install.packages("pheatmap")

#~~~~~~~~~~~~~~~~Default pheatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
distMatheat <- pheatmap::pheatmap(temp, main = "HIV-1/HIV-2 distance matrix plot")

#~~~~~~~~~~~~~~~Customized heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my_color_palette <- colorRampPalette(c("white", "red", "blue"))(100) ##create a color palette

heatmapData <- temp
pheatmap::pheatmap(
  heatmapData,
  color = my_color_palette,  # Custom color palette
  clustering_distance_rows = "euclidean",  # Clustering distance metric for rows
  clustering_distance_cols = "euclidean",  # Clustering distance metric for columns
  cluster_rows = TRUE,  # Whether to cluster rows
  cluster_cols = TRUE,  # Whether to cluster columns
  show_rownames = TRUE,  # Show row names
  show_colnames = TRUE,  # Show column names
  main = " ",  # Main title of the heatmap
  fontsize_row = 9,  # Font size for row labels
  fontsize_col = 9,  # Font size for column labels
  cellwidth = 15,  # Width of heatmap cells
  cellheight = 15  # Height of heatmap cells
)

#####construct neighbor joining tree using ape

tre <- njs(distMat)
tre <- ladderize(tre)
plot(tre, main = "Phylogenenetic analysis of HIV-1 and HIV-2 recombinant sequences")


# Perform hierarchical clustering using a suitable method
clustering <- hclust(distMat)  # Use the appropriate method


#plot the distances
#transform the clustering as a dendrogram object
dendrogram <- as.dendrogram(clustering)

#transform dendrogram into phylo object 
phylotree = as.phylo(clustering)
class(phylotree)

##using plot.phylo 
plot.phylo(phylotree, type = "tidy", main = "Sequence similarity in west African HIV-1 and HIV-2") 

#~~~~~~~~~~~~~~~~~~~~using ggtree~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1<-ggtree(tre, layout = "rectangular") + 
  geom_tippoint(color = "blue", alpha = .5)
p1

p1 + geom_tiplab(hjust = 0, size = 4) + 
  geom_treescale(width = .4) + #define the length of tree
  theme_tree2()

#~~~~~~~~~~~~~~~~~~~~end of script~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
