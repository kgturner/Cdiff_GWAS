#phylogeny for alignment

library(devtools)
install_github("taxize_", "ropensci")

library(taxize)
library(XML)
library(reshape2)

# Define a species list
#from Nucleotide search at https://www.ncbi.nlm.nih.gov/nuccore/
# ("Asteraceae"[Organism] OR ("Asteraceae"[Organism] OR Asteraceae[All Fields])) AND chloroplast[All Fields] AND complete[All Fields] AND genome[All Fields]
# visually scan for length of sequence and statement "complete chloroplast" or "complete plastid"
spplist <- c("Centaurea diffusa","Artemisia frigida","Jacobaea vulgaris",
             "Ageratina adenophora","Lactuca sativa","Helianthus annuus",
             "Guizotia abyssinica","Parthenium argentatum","Helianthus strumosus",
             "Helianthus decapetalus","Helianthus grosseserratus", "Helianthus tuberosus",
             "Helianthus divaricatus", "Helianthus giganteus","Helianthus hirsutus", 
             "Helianthus maximiliani", "Chrysanthemum indicum", "Chrysanthemum x morifolium",
             "Praxelis clematidea")
# Praxelis aka Eupatorium catarium

# 
# Make the tree
# 
# species tree with sub-genera detail (minus three species: parthenium_argentatum, 
# helianthus_sturmosus, chrysanthemum_xmorifolium)
# 
tree <- phylomatic_tree(taxa = spplist, storedtree='smith2011')
tree$tip.label <- taxize_capwords(tree$tip.label)

# # species tree with all of the species 
# tree <- phylomatic_tree(taxa = spplist)
# tree$tip.label <- taxize_capwords(tree$tip.label)


png("Asteraceae_plastome_treeSmith.png", width=600, height = 600, pointsize = 12)

plot(tree, cex=1.1, main="Plastomes in the Asteraceae")
rect(14.8, 15.5, 20.2, 16.4, border="#66CC00",lwd=3)

dev.off()

# legend("bottomleft", legend = c("Weed", "Not a weed"), fill = c("#024DFD", "#ED6500"))