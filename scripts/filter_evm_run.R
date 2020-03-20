setwd("/Volumes/Seagate Backup Plus Drive/Projects/Vienna/Hi5_new/gene_prediction/Compare_gene_predictions/")

evm_out <- read.csv("run/hi5.evm.out", header = F, sep = "\t", stringsAsFactors = F)
evm_gff <- read.csv("run/hi5.evm.gff", header = F, sep = "\t", stringsAsFactors = F)

model_pos <- grep("# EVM", evm_out[,1])
model_end <- c(model_pos[2:length(model_pos)]-1, nrow(evm_out))

model_evidences <- matrix(NA, ncol = 2)

for(mod in 1:length(model_pos)){
  
  c_model <- evm_out[model_pos[mod]:model_end[mod],]
  c_model <- c_model[grep("# EVM", c_model[,1], invert = T),]
  mod_el <- unlist(lapply(strsplit(c_model[,6], ","), length))
  mod_counts <- matrix(data = NA, ncol = 2, nrow = length(mod_el))
  mod_counts[,2] <- mod_el
  mod_counts[,1] <- c_model[,3]
  mod_counts <- rbind(c(paste("Model", mod, sep = "_"), "Nr. of evidences"), mod_counts)
  model_evidences <- rbind(model_evidences, mod_counts)
  
}

model_evidences <- model_evidences[2:nrow(model_evidences),]

model_pos <- grep("Model", model_evidences[,1])
model_end <- c(model_pos[2:18855]-1, nrow(model_evidences))

model_means <- cbind(model_evidences[model_pos,1], NA)

for(mod in 1:length(model_pos)){
  
  model_means[mod,2] <- mean(as.numeric(model_evidences[(model_pos[mod]+1):model_end[mod],2]))

}

bad_models <- grep(F, model_means[,2] > 2)

model_pos <- grep("# EVM", evm_out[,1])
model_end <- c(model_pos[2:18855]-1, nrow(evm_out))

worst_models <- NA

for(bm in bad_models){
  
  c_bad_models <- evm_out[model_pos[bm]:model_end[bm],6]
  mod_el <- unlist(lapply(strsplit(c_bad_models, ","), length))
  if(length(grep(3, mod_el)) == 0){worst_models <- c(worst_models, bm)}
  
}

worst_models <- worst_models[2:length(worst_models)]
most_worst_models <- NA

for(bm in worst_models){
  
  c_bad_models <- evm_out[model_pos[bm]:model_end[bm],6]
  if(length(grep("ev_type", mod_el)) > 0){most_worst_models <- c(worst_models, bm)}
  
}

gene_pos <- grep("gene", evm_gff[,3])
gene_end <- c(gene_pos[2:18855]-1, nrow(evm_gff))

final_genes <- evm_gff[1,]
good_genes <- seq(1,18855,1)
good_genes <- good_genes[is.na(match(good_genes, worst_models))]

good_gff <- evm_gff[1,]

for(g in good_genes){
  
  good_gff <- rbind(good_gff, evm_gff[gene_pos[g]:gene_end[g],]) 
  
}

good_gff <- good_gff[2:nrow(good_gff),]

write.table(good_gff, file = "Hi5.final.genes.gff", sep = "\t", quote = F, col.names = F, row.names = F)



