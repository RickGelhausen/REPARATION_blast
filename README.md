# REPARATION
REPARATION (Ribosome Profiling Assisted (Re-)Annotation of Bacterial genomes) is a pipeline that uses ribosome profiling data for a de novo open reading frame delineation in prokaryotic (bacterial) genomes.

The software was created at VIB-UGent Center for Medical Biotechnology and Lab of Bioinformatics and Computational Genomics (Biobix), University of Gent, belgium.

The REPARATION workflow consist of:

1) Genome wide ribosome profiling profile based on choosen occupancy position.

2) Traverse the genome and generate all possible open reading frames.
	- Generate traning set for random forest model.
	
3) Train a random forest to learn and use it for ORF prediction.

4) Start Site selection.


NB: The ribosome profiling data should be aligned onto the same genome used for the analyzed!

# Requirements
The REPARATION software is primarily written in perl, python2.7 and R. For each dataset the complete REPARATION workflow (from the sam files to final results) took under 45 minutes.
It is also dependant on a set of tool binaries which should all be installed on your system before the pipeline can execute all of its commands.

	Perl packages:
	--------------
		- Perl + v.5.10
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
The tool does not rerquire any special installation outside the prerequisites. 

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
	-g:     Genome fasta file
	-sdir:  The "scripts" directory (avialable within the REPARATION directory), defaults to directory of reparation.pl script
	-db:    fasta database of curated bacteria protein sequences


optional input variables

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

_Ribo-seq_Sense_"psite".bedgraph      Sense bedgraph files for genome wide ribosome profile visualization (psite = 1, 3 or 5)
_Ribo-seq_AntiSense_"psite".bedgraph  Antisense bedgraph files for genome wide ribosome profile visualization (psite = 1, 3 or 5)
_Predicted_ORFs.txt                   List of translated open reading frames predicted by REPARATION
_Predicted_ORFs.bed                   bed file of REPARATION predicted open reading frames
_predicted_ORFs.fasta                 fasta file of predicted translated open reading frame
_plastid_image.png                    Image showing plstid predicted P sites (optional)



# Data

The data sets used for the project can be downloaded from the project website xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



