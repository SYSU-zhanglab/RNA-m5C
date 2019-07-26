Update 2019-7-26:
m5C_caller_multiple.py -- fix a bug in using overall conversion rate

Update 2019-5-19:
m5C_pileup_formatter.py -- fix version error

Step-by-step m5C site calling pipeline (v1.1)
======================================================================================
This is a highly modularized pipeline used in our mRNA m5C sequencing project.


Tested environment
======================================================================================
+--------------------+--------------------------------------+
|System              |CentOS Linux release 7.4.1708 (Core)  |
+--------------------+--------------------------------------+
|CPU                 |8 cores or more                       |
+--------------------+--------------------------------------+
|RAM                 |70 GB or more                         |
+--------------------+--------------------------------------+
|Disk                |50 GB for human metadata              |
|                    +--------------------------------------+
|                    |200 GB for a 30 million PE150 data    |
|                    +--------------------------------------+
|                    |500 GB for a 200 million PE125 data   |
+--------------------+--------------------------------------+
* Windows is not compatible since some modules cannot be installed.


Python version and modules
======================================================================================
Python 2.7.14/2.7.16 with numpy (1.13.3), scipy (0.19.1), pysam (0.12.0.1), Biopython (1.70).
Higher versions of these modules should be compatible.


Tested software and open-source scripts
======================================================================================
+------------------------------+-----------------------------------------------+
|JAVA                          |JRE 1.8.0_131                                  |
+------------------------------+-----------------------------------------------+
|Quality control and           | Cutadapt 1.14                                 |
|formatting                    +-----------------------------------------------+
|                              | Trimmomatic 0.36                              |
+------------------------------+-----------------------------------------------+
|Mapping                       | HISAT2 2.1.0                                  |
|                              +-----------------------------------------------+
|                              | Bowtie 2.2.9/2.3.4.2                          |
|                              +-----------------------------------------------+
|                              | meRanGh/meRanGs (compatible in theory)        |
+------------------------------+-----------------------------------------------+
|Bam processing                |Samtools 1.6                                   |
+------------------------------+-----------------------------------------------+


Customized scripts
======================================================================================
+----------------------------------------+-------------------------------------+
|Name                                    |Usage                                |
+========================================+=====================================+
|Metadata generation                                                           |
+----------------------------------------+-------------------------------------+
|gtf2anno.py                             |Transfer GTF to UCSC annotation      |
|                                        |table (Ensembl format only)          |
+----------------------------------------+-------------------------------------+
|gtf2genelist.py                         |Extract gene/isoform information from|
|                                        |GTF (Ensembl format only)            |
+----------------------------------------+-------------------------------------+
|anno_to_base.py                         |Annotate each base in GTF            |
+----------------------------------------+-------------------------------------+
|anno_to_base_remove_redundance.py       |Remove the redundance                |
+----------------------------------------+-------------------------------------+
|fasta_c2t.py                            |Convert Cs in fasta to Ts            |
+----------------------------------------+-------------------------------------+
|BS_hisat2_index.py                      |Build HISAT2 indexes                 |
+----------------------------------------+-------------------------------------+
|ref_sizes.py                            |Extract the lengths of the references|
+----------------------------------------+-------------------------------------+
|Alignment and pileup                                                          |
+----------------------------------------+-------------------------------------+
|BS_hisat2.py                            |Map reads to the genome              |
+----------------------------------------+-------------------------------------+
|BS_bowtie2.py                           |Map reads to the transcriptome       |
+----------------------------------------+-------------------------------------+
|Bam_transcriptome_to_genome_v1.0.py     |Convert transcriptome alignment to   |
|                                        |genome alignment (BAM to BAM)        |
+----------------------------------------+-------------------------------------+
|concat_bam.py                           |Merge BAM files                      |
+----------------------------------------+-------------------------------------+
|pileup_genome_multiprocessing_v1.4.py   |Split the reference and pileup in a  |
|                                        |multiprocessing manner               |
+----------------------------------------+-------------------------------------+
|m5C_pileup_formatter.py                 |Format pileups                       |
+----------------------------------------+-------------------------------------+
|Call sites                                                                    |
+----------------------------------------+-------------------------------------+
|m5C_caller.py                           |Call m5C sites from the formatted    |
|                                        |pileup files                         |
+----------------------------------------+-------------------------------------+
|m5C_caller_multiple.py                  |Call m5C sites from multiple samples |
+----------------------------------------+-------------------------------------+
|m5C_intersection_single_r1.py           |Identify m5C sites in each sample    |
+----------------------------------------+-------------------------------------+
|m5C_intersection_multi_r1.py            |Identify m5C sites in replicates     |
+----------------------------------------+-------------------------------------+


Installation
======================================================================================
Most of the scripts can run stand-alone. Make sure python modules are installed. You can use ``pip`` to install all the modules required. Python 3 is not compatible recently, however, one can try using python's ``2to3`` util.


Running the pipeline
======================================================================================
Please read ``m5C-BS-seq-step-by-step-computation-protocol-v1.1.pdf`` for details.

If you are using task managers like ``SGE`` or ``SJM``, you can use the script ``m5C_pipeline_generator_qsub.py`` or ``m5C_pipeline_generator_SJM.py`` to generate .sh or .sjm files (for ``SJM``, you need to install ``sjm_tools`` <https://github.com/sysuliujh/Bioinfo-toolkit/tree/master/sjm_tools> first).


Contact
======================================================================================
Please contact ``Jianheng Liu (liujh26@mail2.sysu.edu.cn)`` for questions and bug report.


Citation
======================================================================================
Please cite Huang, T., Chen, W., Liu, J., Gu, N. & Zhang, R. Genome-wide identification of mRNA 5-methylcytosine in mammals. Nature structural & molecular biology, doi:10.1038/s41594-019-0218-x (2019) (https://www.nature.com/articles/s41594-019-0218-x).

