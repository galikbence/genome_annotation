# GENOME ANNOTATION WORKFLOW
Pipeline for eukaryotic genome annotation based on external evidences using [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) and [GeneMark-ES](http://exon.gatech.edu/GeneMark/) gene prediction tools.
#

  ## 1. Collecting external evidences
  
   In order to maximize the genome annotation efficiency you should collect such data that can support a gene model. For example __protein__ and __transcript__ sequences from closely related species. It is highly recommended to download these datasets from reliable source!!! The best if you have __RNA-seq__ data directly from the species of interest. You can assemble the transcripts using the [__Trinity__](https://github.com/trinityrnaseq/trinityrnaseq/wiki) transcriptome assembly tool (either applying the genome-guided method for better results) and annotate the transcripts running the [__Trinotate__](https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinotate-Functional-Annotation) pipeline.
   
   In the next steps we will use the collected/generated data to build gene models for better gene prediction. 
   
   In this workflow, we mainly focus on how to integrate data generated with short read sequencing platforms. Also, you can include evidences, specially whole length transcripts, that were generated using long read sequencing technologies (PacBio, MinION). 

  ## 2. Preparing gene models
  
  In this section we will train AUGUSTUS for another species by the following [tutorial](https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html). The main steps are:
  
     1. COMPILE A SET OF TRAINING AND TEST GENES
     
     2. CREATE A META PARAMETERS FILE FOR YOUR SPECIES
      
     3. MAKE AN INITIAL TRAINING
     
     4. RUN THE SCRIPT optimize_augustus.pl
     
     We will use spliced alignments of de novo assembled transcriptome short reads (RNA-Seq) and 
     spliced alignments of protein sequences against the assembled genomic sequence.
     
     You can find more detailed information in the training manual.
  
  Also, we should prepare hints that the gene prediction tool can incorporate. It  will change the likelihood of gene structures candidates. Therefore, the algorithm will tend to predict gene structures that are in agreement with the hints. In the [readme about AUGUSTUS in the RGASP assessment](http://bioinf.uni-greifswald.de/augustus/binaries/readme.rnaseq.html) a detailed method is described that would produce the necessary inputs for this step. You can find more detailed information reading the AUGUSTUS [tutorial](https://fossies.org/linux/augustus/docs/tutorial/prediction.html#prephints).
  
      Creating hints from RNA-Seq data
      
      Massive amounts of short transcriptome reads first need to be aligned to the genome.  We will assume that we have
      already aligned the reads to the genome and that we have WIG and GFF files.
      
      1. The file coverage.wig contains a coverage graph, that contains for each base in the genome and the number of reads
         alignments that cover the position.
      
      2. The file hints.rnaseq.intron.gff contains likely intron positions, inferred from gaps in the query of the read
         alignments. Together with the intron boundaries the multiplicity (mult) is reported, which counts the number of
         alignments that support the given intron candidate, if there is more than one.
         
      Generate hints about exonic regions from the coverage graph (wig file):
      
        cat coverage.wig | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 \
        --src=W --type=ep --radius=4.5 > hints.rnaseq.ep.gff
        
      Concatenate all hints from all sources into one file:
     
        cat hints.est.gff hints.rnaseq.intron.gff hints.rnaseq.ep.gff > hints.gff
       
We can prepare various gene models and hints file for our genome. We will use these files in Section 5.

  ## 3. Repeat modelling & masking
  
  Repeat modelling and masking repeats is a crucial step in the workflow. We can increase the gene prediction tool speed and efficiency by masking out regions that contains non protein coding elemetns. One of the most popular tools is [__RepeatModeler__](https://github.com/Dfam-consortium/RepeatModeler) (including RepeatMasker).

  ## 4. Predicting ncRNAs & tRNAs
  
  ## 5. Predicting protein coding genes

   ### GeneMark-ES

   ### AUGUSTUS

  ## 6. Combining gene models
  
  ## 7. Annotating protein coding genes
  
    ### BLAST
    
    ### InterProScan (GO, KEGG, domains)
    
    


