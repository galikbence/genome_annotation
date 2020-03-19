# GENOME ANNOTATION WORKFLOW
Pipeline for eukaryotic genome annotation based on external evidences using [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) and [GeneMark-ES](http://exon.gatech.edu/GeneMark/) gene prediction tools.
#

  ## 1. Collecting external evidences
  
   In order to maximize the genome annotation efficiency you should collect such data that can support a gene model. For example __protein__ and __transcript__ sequences from closely related species. It is highly recommended to download these datasets from reliable resource!!! The best if you have __RNA-seq__ data directly from the species of interest. You can assemble the transcripts using the [__Trinity__](https://github.com/trinityrnaseq/trinityrnaseq/wiki) transcriptome assembly tool (either applying the genome-guided method for better results) and annotate the transcripts running the [__Trinotate__](https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinotate-Functional-Annotation) pipeline.
   
   In the next steps we will use the collected/generated data to build gene models for better gene prediction. 
   
   In this workflow, we mainly focus on how to integrate data generated with short read sequencing platforms. Also, you can include evidences, specially whole length transcripts, that were generated using long read sequencing technologies (PacBio, MinION). 

  ## 2. Preparing gene models
  
  In this section we will train AUGUSTUS for another species by the following [tutorial](https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html). 
  
     1. COMPILE A SET OF TRAINING AND TEST GENES
        
        We will use the data from the section 1. You can find more detailed information in the training manual.
     
     2. CREATE A META PARAMETERS FILE FOR YOUR SPECIES
  
  

  ## 3. Repeat modelling & masking

  ## 4. Predicting ncRNAs & tRNAs
  
  ## 5. Predicting protein coding genes

    ### GeneMark-ES

    ### AUGUSTUS

  ## 6. Combining gene models
  
  ## 7. Annotating protein coding genes
  
    ### BLAST
    
    ### InterProScan (GO, KEGG, domains)
    
    


