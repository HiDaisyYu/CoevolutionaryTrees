#########################################
#### Coevolutionary Trees Comparison ####
#########################################

#### 1. Read trees in Newick format using R
## you need two phylogenetic trees nwk files to do below analysis
library(ape)
virustree_Newcick <- read.tree(file = "/Users/yuedana/Desktop/VirusTree.nwk")
hosttree_Newcick <- read.tree(file = "/Users/yuedana/Desktop/HostTree.nwk")
virus_phylo <- as.phylo(virustree_Newcick)
host_phylo <- as.phylo(hosttree_Newcick)
rooted_virus_phylo <- root(virus_phylo,virus_phylo$tip.label[16], resolve.root = TRUE)
rooted_host_phylo <- root(host_phylo, host_phylo$tip.label[16], resolve.root = TRUE)
rooted_host_dend <- chronos(rooted_host_phylo)
rooted_virus_dend <- chronos(rooted_virus_phylo)

#### 2. Comparing trees visually
## Two dendrogram trees were conducted by R package dendextend with script below
> library(dendextend)
> virus_dendex <- as.dendrogram(rooted_virus_dend)
> host_dendex <- as.dendrogram(rooted_host_dend)
> tanglegram(virus_dendex, host_dendex, margin_inner = 8, columns_width = c(5,2,5), main_left = "Coronaviruses", main_right = "Hosts", cex_main = 2)
> den12_untang <- untangle_step_rotate_2side(virus_dendex, host_dendex)
> tanglegram(den12_untang[[1]], den12_untang[[2]], margin_inner = 8, columns_width = c(5,2,5), main_left = "Coronaviruses", main_right = "Hosts")

#### 3. Comparing trees statistically
cor_cophenetic(virus_dendex, host_dendex)
cor_bakers_gamma(virus_dendex, host_dendex)
