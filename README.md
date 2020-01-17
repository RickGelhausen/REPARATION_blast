# MY CHANGES TO REPARATION

I replaced the commericial [usearch --ublast tool](https://www.drive5.com/usearch/) by the open-source [blast tool](https://blast.ncbi.nlm.nih.gov/Blast.cgi). The introduction of blast increases the runtime of reparation, but it allows me to add reparation_blast to bioconda for easy use in conda environments. (avoiding licensing issues caused by usearch)

I also changed the script to soly use bam files instead of a combination of sam and bam files.
Further, I upgraded the python script from python2.7 to python3.6+.

## New dependencies:
- blast 2.7.1+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- python 3.6+


# REPARATION
REPARATION (Ribosome Profiling Assisted (Re-)Annotation of Bacterial genomes) is a pipeline that uses ribosome profiling data for a de novo open reading frame delineation in prokaryotic (bacterial) genomes.

The software was created at VIB-UGent Center for Medical Biotechnology and Lab of Bioinformatics and Computational Genomics (Biobix), University of Gent, belgium.

Link to the paper https://doi.org/10.1093/nar/gkx758

The REPARATION workflow consist of:

1) Genome wide ribosome profiling profile based on choosen occupancy position.

2) Traverse the genome and generate all possible open reading frames.
	- Generate traning set for random forest model.
	
3) Train a random forest to learn and use it for ORF prediction.

4) Start Site selection.


NB: The ribosome profiling data should be aligned onto the same genome used for the analyzed!

REPARATION was developed and tested on a Linux system.


# Requirements

The REPARATION software is primarily written in perl, python2.7 and R. For each dataset the complete REPARATION workflow (from the sam files to final results) takes under 90 minutes.
REAPARATION is dependant on a set of tool binaries which should all be installed on your system before the pipeline can execute all of its commands.

	Perl packages:
	--------------
		- Perl + v.5.10+
		- Getopt::Long
		- Cwd	
		- BioPerl
		- POSIX
		
	Python packages:
	--------------
		- plastid (https://plastid.readthedocs.io/en/latest/)
		
	R packages:
	--------------
		- ggplot2
		- randomForest
		- PRROC
		- ROCR
		- SiZer

	Tool binaries:
	--------------
		- prodigal[*] (https://github.com/hyattpd/Prodigal)
		- glimmer3 (https://ccb.jhu.edu/software/glimmer/)
		- USEARCH[*] (http://www.drive5.com/usearch/download.html)
		- R (http://www.r-project.org/)
		- samtools
		
		[*] The usearchX.Y.Z and prodigalX.Y.Z executatble should be renamed to "usearch" and "prodigal" respectively.
	
	The tool binary paths should be included in the $PATH variable or within the script/bin folder of the tool.

	
# Install

REPARATION does not rerquire any special installation requirements. To install REPARATION, clone/download the github repository into a directory of your choice.

Binaries:

If prodigl/glimmer3 is not installed on your system then REPARATION will use the version avialable in the script directory. Under Linux you must ensure that you have read and execute permissions for the binary file. If needed, use the chmod command to set the execute bit, e.g.:

chmod +x $PATH/REPARATION/scripts/bin/prodigal or chmod 755 $PATH/REPARATION/scripts/bin/prodigal

chmod +x $PATH/REPARATION/scripts/bin/glimmer/glimmer3

chmod +x $PATH/REPARATION/scripts/bin/glimmer/build-icm



# Galaxy Installation

1) Add tool definition files in the galaxy tools directory
2) add tool section in the tool_conf.xml file
	```xml
	<section name="REPARATION" id="reparation">
		<tool file="REPARATION/reparation.xml" />
	</section>
	```
# Usage

Usage: ./reparation.pl -sam riboseq_alignment_files_sam_format -g genome_fasta_file -sdir REPARATION_scripts_directory -db curated_protein_db(fasta) [options]


Mandatory input variables

	-sam:   Ribosome alignment file (sam)
	-g:     Genome fasta file. This should be the same genome fasta file used in the alignment of the Ribo-seq reads.
	-sdir:  The "scripts" directory (avialable within the REPARATION directory), defaults to directory of reparation.pl script
	-db:    fasta database of curated bacteria protein sequences


Optional input variables

	-gtf:   GTF genome annotation file
	-wdir:  working directory (defualts to current directory)
	-en:    Experiment name
	-p:     Ribosome profiling read p site assignment strategy, 1 = plastid P-site estimation ((default), 3 = 3 prime of read, 5 prime of the read
	-mn:    All ribosome profiling reads shorter than these values are eliminated from the ananlysis (default = 22)
	-mx:    All ribosome profiling reads longerer than these values are eliminated from the ananlysis (default = 40)
	-mo:    Minimum length of open reading frame considered for prediction (default = 30 value in nucleotides)
	-mr:    Open reading frames with less than these number of ribosome profiling reads are eliminated from analysis (default = 3)
	-ost:   Start region length in nucleotides (default = 45nts). This value is used to calculate features specific to the start region.
	-osp:   Stop region length in nucleotides (default = 21nts). This value is used to calculate features specific to the stop region.
	-osd:   Distance of Shine dalgarno sequence from start codon (defualt = 5nts). 
	-seed:  Shine dalgarno sequence (default = GGAGG). The shine dalgarno sequence used for Ribosome binind site energy calculation.
	-sd:    Use ribosome binding site (RBS) energy in the open reading frame prediction (Y = use RBS energy (default), N = donot use RBS engergy)
	-id:    Minimum identify score for BLAST protein sequence search (defualt = 0.75)
	-ev:    maximum e-vlaue for BLAST protein sequence search (default = 1e-5)
	-pg:    program for initial positive set generation (1 = prodigal (default), 2 = glimmer)
	-cdn:   Comma separted list of start codons (default = ATG,GTG,TTG)
	-ncdn:  Start codon for negative set (default = CTG)
	-pcdn:  Start codon for positive set (default = ATG,GTG,TTG). Should be a subset of the standard genetic code for bacteria
    -gcode: Genetic code to use for initail positive set generation. Valid when -pg is 1. (default = 11, takes value between 1-25)
    -by:    Flag to determine if prodigal should bypass Shine-Dalgarno trainer and force a full motif scan (default = N). (Y = yes, N = no) Valid only for -pg 1
    -score: Random forest classification probability score threshold to define as ORF are protein coding, the minimum  (defualt is 0.5)


Output files

_Ribo-seq_Sense_"psite".bedgraph      	Sense bedgraph files for genome wide ribosome profile visualization (psite = 1, 3 or 5)

_Ribo-seq_AntiSense_"psite".bedgraph  	Antisense bedgraph files for genome wide ribosome profile visualization (psite = 1, 3 or 5)

_Predicted_ORFs.txt                   	List of translated open reading frames predicted by REPARATION

_Predicted_ORFs.bed                   	bed file of REPARATION predicted open reading frames

_predicted_ORFs.fasta                 	fasta file of predicted translated open reading frame

_plastid_image.png			Image showing plstid predicted P sites (optional)

_PR_ROC_curve.pdf                       Precision-Recall and ROC curve plots

_metagene_profile.pdf                   Metagene profiles around the start and stop of ORFs in positive set

_Scurve.pdf                             Sigmoid curve with estimated thresholds

_variable_importance.pdf                variable importance plot

_PR_ROC_curve.pdf		       	Plots for the Precision-Recall and ROC curve to evalaute model performance

_metagene_profile.pdf		       	Metagene profile around the start and stop of ORFs in the positive set

_Scurve.pdf			       	Plot of the Sigmoid curve showing the estimated minimum thresholds

_variable_importance.pdf	       	Varible importance (Gini) of the faetures used in the model


# Data

The data sets used for the project can be downloaded from http://www.biobix.be/reparation/data/

# More information

For more information about REPARATION contact: elvis.ndah@gmail.com or Gerben.Menschaert@UGent.be
