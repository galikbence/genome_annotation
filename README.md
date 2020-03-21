# GENOME ANNOTATION WORKFLOW
Workflow for eukaryotic genome annotation based on external evidences using [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) and [GeneMark-ES](http://exon.gatech.edu/GeneMark/) gene prediction tools.
#

   ![workflow](https://user-images.githubusercontent.com/47104867/77142207-ac294c80-6a7f-11ea-8de6-6667d8919472.png)

__Figure 1.__ Genome annotation workflow.
#

The workflow contains the following steps:

  1. Collecting external evidences
  2. Preparing gene models
  3. Repeat modelling & masking
  4. Predicting tRNAs
  5. Predicting protein coding genes
  6. Combining gene models
  7. Annotating protein coding genes
  
 #
 
  ## 1. Collecting external evidences
  
  We assume you already sequenced and assembled your genome.
  
   In order to maximize the genome annotation efficiency you should collect such data that can support a gene model. For example __protein__ and __transcript__ sequences from closely related species. It is highly recommended to download these datasets from reliable source!!! The best if you have __RNA-seq__ data directly from the species of interest. You can assemble the transcripts using the [__Trinity__](https://github.com/trinityrnaseq/trinityrnaseq/wiki) transcriptome assembly tool (either applying the genome-guided method for better results) and annotate the transcripts running the [__Trinotate__](https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinotate-Functional-Annotation) pipeline.
   
   In the next steps we will use the collected/generated data to build gene models for more accurate gene prediction. 
   
   In this workflow, we mainly focus on how to integrate data generated with short read sequencing platforms. Also, you can include evidences, specially whole length transcripts, that were generated using long read sequencing technologies (PacBio, MinION). 

  ## 2. Preparing gene models and other evidences
  
  In this section we will train AUGUSTUS for another species by the following [tutorial](https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html). Also we will prepare other evidences.
  
  __Retraining AUGUSTUS__
  
  The main steps are:
  
     1. COMPILE A SET OF TRAINING AND TEST GENES
     
     2. CREATE A META PARAMETERS FILE FOR YOUR SPECIES
      
     3. MAKE AN INITIAL TRAINING
     
     4. RUN THE SCRIPT optimize_augustus.pl
     
     We will use spliced alignments of de novo assembled transcriptome short reads (RNA-Seq) and 
     spliced alignments of protein sequences against the assembled genomic sequence.
     
     You can find more detailed information in the training manual.
  
  Also, we should prepare hints that the gene prediction tool can incorporate. It  will change the likelihood of gene structures candidates. Therefore, the algorithm will tend to predict gene structures that are in agreement with the hints. In the [readme about AUGUSTUS in the RGASP assessment](http://bioinf.uni-greifswald.de/augustus/binaries/readme.rnaseq.html) a detailed method is described that would produce the necessary inputs for this step. You can find more detailed information reading the AUGUSTUS [tutorial](https://fossies.org/linux/augustus/docs/tutorial/prediction.html#prephints).
  
 __Creating hints from RNA-Seq data__
      
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
        
        
 __Creating hints from cDNA__
 
        blat -minIdentity=92 <genome.fa> <cdna.fa> <cdna.psl>
   
        pslCDnaFilter -maxAligns=1 cdna.psl cdna.f.psl

        blat2hints.pl --in=cdna.f.psl --out=hints.cdna.gff
        
      
 __Concatenate all hints from all type sources into one file:__
     
        cat hints.cnda.gff hints.rnaseq.intron.gff hints.rnaseq.ep.gff > hints.gff
        
      
       
We can prepare various gene models and hints file for our genome. We will use these files in Section 5.

  ## 3. Repeat modelling & masking
  
  Repeat modelling and masking repeats is a crucial step in the workflow. We can increase the gene prediction tool speed and efficiency by masking out regions that contains non protein coding elements. One of the most popular tools is [__RepeatModeler__](https://github.com/Dfam-consortium/RepeatModeler) (including RepeatMasker). 
  
  The main steps are:
  
    1. Create a Database for RepeatModeler

       <RepeatModelerPath>/BuildDatabase -name genome_of_interest genome_of_interest.fasta

    2. Run RepeatModeler

       <RepeatModelerPath>/RepeatModeler -database genome_of_interest -pa 20 -LTRStruct >& run.out &
       
    3. Interperting the results
    
       At the succesful completion of a run, two files are generated:
       
         <database_name>-families.fa  : Consensus sequences
         <database_name>-families.stk : Seed alignments
         
    4. Making repeat library
    
       <RepeatMaskerPath>/RepeatMasker -lib <database_name>-families.fa mySequence.fa
       
    5. Predicting repeats
    
       <RepeatMaskerPath>/RepeatMasker -small -gff --species <query species> --lib [filename for custom library] yourgenome.fasta
  
  A few gene prediction tool can recognize masked regions if your genom sequence is upper case and the masked regios are in lower case. Therefore we highly recommend to use `-small` option. Later, we will need the repeats positions information in gff fromat (use option `-gff`) in Section 6. 
  
  Eukaryotic genomes can be huge and repeat masking can take a lot of time (also memory). If the `-pa(rallel)` option is not working you can speed up this step by splitting up your genome into separate FASTA files and you can mask these files parallel. At the end you can simply combine the masked FASTA and GFF files into one.

  ## 4. Predicting tRNAs
  
  There are not only protein coding genes in your _genome of interest_. If you want to see the whole picture you should predict __tRNA__ coding genes too. [tRNAscan-SE](https://github.com/biopro/genix/tree/master/bin/tRNAscan-SE) is one of the best tools for tRNA gene predictions. It is easy to run, fast and it has a [webserver](http://lowelab.ucsc.edu/tRNAscan-SE/).
  
  Simple usage (see more options in the [manual](https://github.com/biopro/genix/blob/master/bin/tRNAscan-SE/MANUAL)):
  
      tRNAscan-SE [-options] <FASTA file(s)>
   
  ## 5. Predicting protein coding genes
  
  This is the core part of the workflow. We will use an **evidence driven _ab initio_ gene prediction approach** running several gene predictor tools. Since we want to predict thousand of protein coding genes we automatize the process but we don't have time to refine all gene models by hand. Therefore, we will use more than one tool and in the next section we will combine these results. Also, each algorthim has it's advandage. 
  
  We will run the following algorithms:

   ### GeneMark-ES
   
 New genomes can be analyzed by [GeneMark-ES](http://exon.gatech.edu/GeneMark/) applying un/supervised self-training which is an important feature of the algorithm. Also, it can take external evidences (RNA or protein) mapped to genome. So, we can use the previously prepared files (see Section 2.). We can run the `gmes_petap.pl` perl script to do the trainig and prediction steps in one. See the [documentation](https://wiki.gacrc.uga.edu/wiki/GeneMarkES-Teaching) for more information.
   
         gmes_petap.pl --soft_mask --ES --evidence hints.gff --cores <number of cores> --sequence genome_of_interest.fasta
         
  The ouput is a GTF file. However, the program can predict incomplete genes but we are not interesed in these gene models. We can filter out using the `filter_genemark.R` script that you can find in the [repository](https://github.com/galikbence/genome_annotation/tree/master/scripts). Next we should convert GTF into GFF using [gffread](https://github.com/gpertea/gffread).
  
  You should run GeneMark-ES only once!

   ### AUGUSTUS

AUGUSTUS is a program that predicts genes in eukaryotic genomic sequences. It can be run on web server or run locally. It has 2 mandatory arguments. The query file and the species. You can find more details about how to run the tool in the [manual](https://github.com/Gaius-Augustus/Augustus/blob/master/docs/RUNNING-AUGUSTUS.md).

You can run AUGUSTUS several times using the basic gene model, the gene model of the closest related species and the gene model with hints that we prepared in Section 2.

Exmple runs (each case we will predict only complete genes on both strands):

      #Run 1
       
       augustus --strand=both --genemodel=complete --species=[related_speices] --gff3=on --codingseq=on --outfile=[out_file] genome_of_interest.fasta
       
      #Run 2
       
       augustus --strand=both --genemodel=complete --species=[related_speices] --gff3=on --codingseq=on --hintsfile=hints.gff --extrinsicCfgFile=extrinsic.ME.cfg --outfile=[out_file] genome_of_interest.fasta
       
      #Run 3
       
       augustus --strand=both --genemodel=complete --species=generic --gff3=on --codingseq=on --hintsfile=hints.E.gff --extrinsicCfgFile=extrinsic.ME.cfg --outfile=[out_file] genome_of_interest.fasta
       
At the end of these analyses we will use the `GFF` files from each gene prediction run (GeneMark-ES and AUGUSTUS). You may have to uniform the `GFF` file format amongst the different runs.

  ## 6. Combining gene models
  
In this section we will use all external evidences, repeats and gene prediction results in order to merge the different gene models into one.

The [EVidenceModeler](https://evidencemodeler.github.io) (EVM) software combines gene predictions and protein and/or transcript alignments (as external evidences) into __weighted__ consensus gene models. EVM provides a flexible and intuitive framework for combining diverse evidence types into a single automated gene structure annotation system.

Inputs to EVM include: 
   - the genome sequence 
   - gene predictions in GFF3 format
   - alignment data in GFF3 format 
   - list of numeric weight values to be applied to each type of evidence
   
The weights can be configured manually!

   ### Preparing the inputs

   #### 1. Gene predictions
   
   We should merge all `GFF` files from all gene predictions into one file. 
   
      #Example
      
       cat prediction_1.gff prediction_2.gff prediction_3.gff prediction_N.gff > gene_predictions.gff
      
   If the separate `GFF` files contains headers (lines statring whit "#") you sould delete them. Other important thing is the __source__ column (the 2nd one) of `GFF` files should be __unique__ corresponding to each run.
   
       #Example GFF files
       
        Scaffold2	AUGUSTUS_1	gene	7345	7746	1	-	.	ID=g1
        Scaffold2	AUGUSTUS_2	gene	7345	7746	1	-	.	ID=g1
        Scaffold2	AUGUSTUS_3	gene	7345	13968	0.07	-	.	ID=g1
        Scaffold2	GeneMark.hmm	gene	7345	7746	1	-	.	ID=g1
   
   #### 2. Alingments for evidences
   
   We can prepare cDNA (assembled transcripts from RNA-seq data) alignments using [GMAP/GSNAP](https://github.com/juliangehring/GMAP-GSNAP)
   
        #Creating genome index
         gmap_build -d <genome> [-k <kmer size>] <fasta_files...>
        
        #Mapping transcripts to the genome
         gmap -d <genome> -B 5 -t 10 -f 3 transcripts.fasta
         
 Rename the output file as `cdna.gff`!
         
   Also, we can use protein based evidences generated with [Exonerate](https://github.com/nathanweeks/exonerate) and [Scipio](https://www.webscipio.org).
   
         #Exonerate example
          exonerate --model protein2genome query_proteins.fasta target_genome.fasta
          
         #Scipio example
          scipio.pl [<options>] <target genome> <query protein>

See more option in the manuals ([Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-user-guide), [Scipio](https://www.webscipio.org/help/scipio)).

Merge then rename the output file as `proteins.gff`! Again, the __source__ column (the 2nd one) of `GFF` files should be __unique__ corresponding to each run. 

   #### 3. Repeats
   
We can use the `GFF` output from repeat masing (see Section 3.) Rename the output as `repeats.gff`.

   #### 4. Weights
   
We can easily creat the weights file which has three columns: 
   - evidence class 
   - type, correspoding to the name of the tool/run that was used to generate it
   - weight
   
The class parameter can be one of the following: 
   - ABINITIO_PREDICTION 
   - PROTEIN
   - TRANSCRIPT
   
These are the only input types accepted by EVM currently. An example weight file looks like:

         TRANSCRIPT	GMAP	10
         PROTEIN	exonerate	5
         PROTEIN	Scipio	5
         ABINITIO_PREDICTION	AUGUSTUS_1	1
         ABINITIO_PREDICTION	AUGUSTUS_2	1
         ABINITIO_PREDICTION	AUGUSTUS_3	1
         ABINITIO_PREDICTION	GeneMark.hmm	1
   
   ### Running EVM
   
We have all the necessary inputs for running EVidenceModeler. 

The main steps are:

- Partitioning the Inputs

      perl $EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome sotfmasked_genome.fasta \
           --gene_predictions gene_predictions.gff --protein_alignments proteins.gff \
           --transcript_alignments cdna.gff \
           --repeats repeats.gff \
           --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out

 - Generating the EVM Command Set
      
       perl ${EVM_HOME}/EvmUtils/write_EVM_commands.pl --genome sotfmasked_genome.fasta --weights weights.txt \
            --gene_predictions gene_predictions.gff --protein_alignments proteins.gff \
            --transcript_alignments cdna.gff \
            --repeats repeats.gff \
            --output_file_name evm.out  --partitions partitions_list.out >  commands.list

       perl $EVM_HOME/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log

- Combining the Partitions

      perl $EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

- Convert the results to GFF3 Format

      perl $EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome sotfmasked_genome.fasta
      
- Merge GFF files
  
      cat <Name of the contigs/scaffolds>/evm.out > evm.all.out
      cat <Name of the contigs/scaffolds>/evm.out.gff3 > evm.out.all.gff

Or you can use the [runEVidenceModeler.sh](https://github.com/galikbence/genome_annotation/blob/master/scripts/runEvidenceModeler.sh) script. Before running it you may have change it a little bit!

At the end of this process we can filter the results based on the number of evidences that support a gene model. For example we get rid of those models that have only 1 _ab initio_ evidence. You can use the [filter_evm_run.R]() script with basic filtering options. The input are the `evm.out.all.gff` and `evm.all.out` files.
  
Finally we will extract the CDS and protein sequences that we will use in the annotation step.

In this workflow we don't care about the alternative splicing events. However, it can add more to the annotation.

  ## 7. Annotating protein coding genes
  
 Now we have all the protein coding genes from Section 6. We will annotate structural and functional features of the coding sequences (or you can use the proteins). One of the most popular approach is to annotate the gene of interes based on similarity unsing [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi). If you have to  align thousands of sequnces you can use the local version and you can prepare a custom database by colleting all __well annotated__ proteins/transcripts from closely related species or you can use the [UniProtKB/Swiss-Prot](https://www.uniprot.org/statistics/Swiss-Prot) database. You can align DNA to DNA or PROTEIN to DNA.
 
 See more in the BLAST [documentation](https://www.ncbi.nlm.nih.gov/books/NBK279680/) how to use the command line version.
 
 However, BALST has good balance of sensitivity and speed, it is flexible and reliable BUT the quality of the results depends on the database. It can happen there are a lot of unannotated and hypothetical elements in your databese. Also, if your databes is good you will have a lot of elemets that have poor matches etc... 
 
 Therefore, we highly recommend to scan your sequence for matches against the protein signature databases. For this porupse [InterPro](https://www.ebi.ac.uk/interpro/) provides functional analysis of proteins by classifying them into families and predicting domains and important sites. You can do quick [search](https://www.ebi.ac.uk/interpro/search/sequence/) with limited number of amino acids (40,000). If you want to analyze thousands of sequences you can download the [InterProScan](https://www.ebi.ac.uk/interpro/download/) tool. It has a github [repository](https://github.com/ebi-pf-team/interproscan) but for more information on downloading, installing and running it please see the [wiki](https://github.com/ebi-pf-team/interproscan/wiki) page.
  
   ### BLAST
    
   ### InterProScan (GO, KEGG, domains)
    
    


