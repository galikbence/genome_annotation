#!/usr/local/bin/R

#Load packages
library("argparser")

#Options
parser <- arg_parser(description='Process commandline arguments')
parser <- add_argument(parser, arg=c("--gtf"), 
                       help = c("GeneMark-ES results in GTF format"),
                       type = c("character"))


args = parse_args(parser)

#Load GeneMark-ES result
genemark_gtf <- read.csv(args$gtf, sep = "\t", stringsAsFactors = F, header = F)

#Selecting genes on the PSITIVE and NEGATIVE strand 
pos_strand <- genemark_gtf[grep("-", genemark_gtf[,7], invert = T),]
neg_strand <- genemark_gtf[grep("-", genemark_gtf[,7]),]

#Get START and STOP codons positions
pos_start_pos <- grep("start_codon", pos_strand[,3])
pos_stop_pos <- grep("stop_codon", pos_strand[,3])
neg_start_pos <-grep("start_codon", neg_strand[,3])
neg_stop_pos <- grep("stop_codon", neg_strand[,3])

#Filter complete genes in "+" strand
pos_genes <- sort(c(pos_start_pos, pos_stop_pos))
pos_ord <- pos_strand[pos_genes,]
good_pos_genes <- pos_ord[grep(TRUE, duplicated(pos_ord[,9])),9]

#Filter complete genes in "-" strand
neg_genes <- sort(c(neg_start_pos, neg_stop_pos))
neg_ord <- neg_strand[neg_genes,]
good_neg_genes <- neg_ord[grep(TRUE, duplicated(neg_ord[,9])),9]

#Combine final list of good genes
good_genes <- c(good_pos_genes, good_neg_genes)
good_gtf <- genemark_gtf[grep(FALSE, is.na(match(genemark_gtf[,9], good_genes))),]

#Write ouput
write.table(good_gtf, file = "good_genemark.gtf", sep = "\t", quote = F, col.names = F, row.names = F)
