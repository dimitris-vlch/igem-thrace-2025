library(msa)
library(seqinr)
library(ape)

hemoSeq <- readAAStringSet("sequences.fasta")

hemoAln <- msa(hemoSeq)
hemoAln2 <- msaConvert(hemoAln, type="seqinr::alignment")
d <- dist.alignment(hemoAln2, "identity")

d = as.matrix(d)
d = d[rowSums(is.na(d)) == 0, colSums(is.na(d)) == 0, drop = FALSE]
d = as.dist(d)
write.table(as.matrix(d),"dd")

hc = hclust( d , method = "complete")
plot(hc, hang = -1, cex = 0.4)


tree <- read.tree("sequences.dnd")
plot(tree, main="Phylogenetic Tree")


tree <- read.tree("sequences2.dnd")
plot(tree, main="Phylogenetic Tree")
