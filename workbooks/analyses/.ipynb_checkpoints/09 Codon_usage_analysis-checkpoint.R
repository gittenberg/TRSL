# script adapted from https://github.com/mariodosreis/tai

#setwd("~/tai-master/misc")
setwd("~/git/TRSL/workbooks/analyses/tAI_misc")

require("tAI")
#yeast.trna = scan("yeast.trna")  # wrong according to https://github.com/smsaladi/tAI/issues/2
yeast.trna = scan("sc_abundance.trna")
yeast.ws = get.ws(tRNA=yeast.trna, sking=0) # superkingdom, 0 indicating Eukaryota
yeast.m = matrix(scan("yeast.m"), ncol=61, byrow=TRUE)
#yeast.m = matrix(scan("ecolik12.m"), ncol=61, byrow=TRUE)  # test (to replicate https://github.com/smsaladi/tAI/issues/2, set sking=1)
yeast.m = yeast.m[,-33]

yeast.tai = get.tai(yeast.m, yeast.ws)
hist(yeast.tai)

# annotate tAI by name

#biocLite("Biostrings")  # to install Biostrings
library(Biostrings)
dna = readDNAStringSet("orf_coding.fasta")

# get names only (from https://stackoverflow.com/questions/20428742/select-first-element-of-nested-list)
s = strsplit(names(dna), "[ \t]+") # split names by tab/space
gene.names = lapply(s, function(l) l[[1]])

# name the list (https://stackoverflow.com/questions/17842705/creating-a-named-list-from-two-vectors-names-values)
names(yeast.tai) = gene.names

library(RJSONIO)
yeast.json = toJSON(list(yeast.tai)[[1]])
write(yeast.json, file="yeast.json")
