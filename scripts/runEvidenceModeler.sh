#!/bin/bash
EVM_HOME=/Path/to/EVidenceModeler

echo "Partitioning the Inputs"
perl $EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome sotfmasked_genome.fasta \
     --gene_predictions gene_predictions.gff --protein_alignments proteins.gff \
     --transcript_alignments cdna.gff \
     --repeats repeats.gff \
     --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out

echo "Generating the EVM Command Set"
perl ${EVM_HOME}/EvmUtils/write_EVM_commands.pl --genome sotfmasked_genome.fasta --weights weights.txt \
      --gene_predictions gene_predictions.gff --protein_alignments proteins.gff \
      --transcript_alignments cdna.gff \
      --repeats repeats.gff \
      --output_file_name evm.out  --partitions partitions_list.out >  commands.list

perl $EVM_HOME/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log

echo "Combining the Partitions"
perl $EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out


echo "Convert to GFF3 Format"
perl $EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome sotfmasked_genome.fasta

echo "Merge GFF files"
cat <Name of the contigs/scaffolds>/evm.out > evm.all.out
cat <Name of the contigs/scaffolds>/evm.out.gff3 > evm.out.all.gff
