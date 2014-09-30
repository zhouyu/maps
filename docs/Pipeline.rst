==============================
Pipeline of MAPS data analysis
==============================
|

Here is a simple list of steps to do MAPS data analysis.

Decoding sequencing data
------------------------
|

We generally multiplex 12 or more samples in one lane of HiSeq2000 Ilumina sequencing run. Each sample has a unique barcode sequence. After getting the Fastq file (f_fastq), we can decode and split the data into multiple Fastq files per sample with script *decode.py*.

  1. Create a barcode configuration file as f_barcode (see an example below)
  2. Get the start position of barcode in sequencing read (37 by default)
  3. Choose maximum number of mismatches (0 by default)
  4. Run script *decode.py* to decode samples (see an example below) 
  5. Run script *fastq2dropend.py* to keep only the target sequence by providing the length to drop for each read

An example of f_barcode: 4 samples with barcode sequence in 1st column, and sample in 2nd column (separated by \\t, tab).

.. code-block::

  ACTG	WT1
  GGAT	KO1
  AGGA	WT2
  CGTC	KO2

Decoding usage example:

.. code-block::

  decode.py [--startpos=37 –mismatch=0 –outprefix lane1] f_fastq f_barcode

This will generate Fastq files named as lane1_WT1.fastq, etc. 

.. code-block::

  fastq2dropend.py -d 7 -o WT1clean.fastq WT1.fastq

This will drop last 7 nt of reads in WT1.fastq and write reads to file WT1clean.fastq.


Mapping sequencing reads
------------------------
|

The user could choose different mappers to do mapping from reads to the locations where the reads are from. We generally use Bowtie to map reads to the genome, for example human genome hg18.  Here are the basic steps.

  1. Build genome index with program *bowtie-build* from sample genome sequence
  2. Remove 3'-adaptor and/or polyA sequences at the end of read with program `cutadapt <https://github.com/marcelm/cutadapt>`_
  3. Keep reads longer than 18 nt with simple script
  4. Map reads to the reference and keep unique hits (see example below)

.. code-block::

  bowtie EBWT --trim5 4 -l25 -n2 -k1 --best --strata

EBWT is the reference index. We skip the first 4 nt of reads in mapping, as the beginning nucleotides from random primers may not match template, and allow two mismatches in the remaining first 25 nt seed sequence.


Create reference genes
----------------------
|

We generally use annotated genes from `UCSC genome browser <http://genome.ucsc.edu/>`_ to compute their expression. As MAPS method theoretically generate exonic reads close to the gene's 3-end, we normally only count reads close to the 3'-end. And also, one gene may have multiple annotated 3'-ends, we count all of them for computing the gene's expression level. Users could compute expression on different 3'-ends separately to investigate APA switch problem. Below are basic steps to create the "landing exons" of reads for reference genes.


  1. Download refGene genomic coordinates in BED format with UCSC table browser
  2. Rename the gene IDs (column 4) to unique names (for example row id)
  3. Cluster isofoms/transcripts of genes with program *clusterGenes*
  4. Create the counting intervals for each gene with script *gene2land.py*

.. code-block::

  gene2land.py --lmrna 300 -o f_land f_cluster f_gene


f_gene is the BED file from step 2 and f_cluster is from step 3. The intervals ("landing exons" of reads) is set by option -lmrna. The output file f_land is in GTF format.


Filtering reads resulting from internal priming
-----------------------------------------------
|

We use a heuristic filter based on previous reported motifs to remove potential internal priming events. The downstream sequence (up to 300 nt) of each read is checked on the presence of polyA stretches (consecutive 8 As or 9 As in a 10 nt window). The polyA stretches can be scanned by the script *scan4astretch.py* (see usage below) across the sequences of all annotated genes or the genome once and saved for later use. The output BED file can be compared with the mapped reads by the *intersectBed* program in `BEDTools <https://github.com/arq5x/bedtools2>`_ to filter out those potential internal priming events.

.. code-block::

  scan4astretch.py -o polyA.bed --ival gene.bed genome_dir

The directory genome_dir contains the Fasta files of the sample genome sequences. The search can be limited in annotated genes with option --ival by provided the BED file of gene intervals. 


Compute gene expression
------------------------
|

Finally, we use the script *land2exp.py* to count reads and calculate gene expression in RPKM. Below is a usage example.

.. code-block::

  land2exp.py f_tagbam f_land -u cluster_id -o f_count


f_tagbam is the file of mapped reads in BAM format. The option -u/--unit sets the computation mode for counting reads for a clustered gene (cluster_id)  or a transcript (transcript_id). 

