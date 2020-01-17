#!/usr/bin/env perl

#####################################
##	REPARARTION: Ribosome Profiling Assisted (Re-) Annotation of Bacterial genomes
##
##	Copyright (C) 2017 Elvis Ndah
##
##	This program is free software: you can redistribute it and/or modify
##	it under the terms of the GNU General Public License as published by
##	the Free Software Foundation, either version 3 of the License, or
##	(at your option) any later version.
##
##	This program is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##	GNU General Public License for more details.
##
##	You should have received a copy of the GNU General Public License
##	along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##	contact: elvis.ndah@gmail.com
#####################################

use strict;
use warnings;
use diagnostics;

use Getopt::Long qw< :config auto_version >;
use Pod::Usage;
use Cwd;
use File::stat;
use File::Basename;

$main::VERSION = "1.0.7";
my $startRun = time();	# track processing time


########################
##	Usage: ./reparation.pl -bam riboseq_alignment_files_bam_format -g bacteria_genome_fasta_file -sdir scripts_directory -db curated_protein_db(fasta)
########################

#----------------------------------------------------
#			VARIBLES
#----------------------------------------------------


my %translationHash =
	(GCA => "A", GCG => "A", GCT => "A", GCC => "A",
     TGC => "C", TGT => "C",
     GAT => "D", GAC => "D",
     GAA => "E", GAG => "E",
     TTT => "F", TTC => "F",
     GGA => "G", GGG => "G", GGC => "G", GGT => "G",
     CAT => "H", CAC => "H",
     ATA => "I", ATT => "I", ATC => "I",
     AAA => "K", AAG => "K",
     CTA => "L", CTG => "L", CTT => "L", CTC => "L", TTA => "L", TTG => "L",
     ATG => "M",
     AAT => "N", AAC => "N",
     CCA => "P", CCT => "P", CCG => "P", CCC => "P",
     CAA => "Q", CAG => "Q",
     CGA => "R", CGG => "R", CGC => "R", CGT => "R",
     AGA => "R", AGG => "R",
     TCA => "S", TCG => "S", TCC => "S", TCT => "S",
     AGC => "S", AGT => "S",
     ACA => "T", ACG => "T", ACC => "T", ACT => "T",
     GTA => "V", GTG => "V", GTC => "V", GTT => "V",
     TGG => "W",
     TAT => "Y", TAC => "Y");


# Mandatory variables
my $genome;				    # Prokaryotic genome file in fasta format
my $bam_file;				# Ribosome profiling alignment file [bam format]
my $blastdb;				# protein blast database in fasta format
my $dirname = dirname(__FILE__);	# get tool directory for default script directory
my $script_dir = $dirname."/scripts";	# Directory where the script files are stored (defaults to the current directory)


# Input variables
my $threads = 1;			# number of threads used by the USEARCH tool
my $workdir;				# working directory to store files (defaults to current directoy)
my $experiment;				# Experiment name
my $gtf;				    # Genome annotation file (gtf) [if avialable]
my $occupancy = 1;			# p-site of on reads (1 = plastid estimated p-site (default), 3 = 3 prime end of read  and 5 = 5 prime end of read)
my $min_read_len = 22;		# Minimum RPF read length
my $max_read_len = 40;		# Maximum RPF read length
my $MINORF = 30;			# Minimum ORF length
my $MINREAD = 3;			# Only ORFs with at least this number of RPF reads
my $OFFSET_START= 45;		# Offset at the start of the ORF
my $OFFSET_STOP	= 21;		# Offset at stop of the ORF
my $OFFSET_SD = 5;			# distance uptream of start codon to start search for SD sequence and position
my $SEED = "GGAGG";			# The seed shine dalgano sequence
my $USESD = "Y";				# Flag to determine if RBS energy is included in the predictions [Y = use RBS, N = do not use RBS]
my $identity = 0.75;		# identity threshold for comparative psotive set selection
my $evalue = 1e-5;			# e value threshold for comparative psotive set selection
my $pgm = 1;				# Program to generate positive set prodigal=1, glimmer=2
my $start_codons = "ATG,GTG,TTG";	# Comma seperated list of possible start codons
my $start_codon_nset = "CTG";		# Start codon for the negative set
my $start_codon_pset;		# Start codon for the positive set. Defualts to start codons set. Should be a subste of the standard genetic code for bacteria
my $genetic_code = 11;		# the genetic code [1-25] that determines the allowed start codons
my $seedBYpass = "N";       # Bypass Shine-Dalgarno trainer and force a full motif scan (default = N(o)). Valid only for -pg 1
my $score = 0.5;           # Random forest classifier threshold to classify ORF as protein copding (defualt is 0.5).
my $output_folder = "reparation"; # Name of the output folder for the results.

# Output files
my $bedgraphS;
my $bedgraphAS;
my $predicted_ORFs;
my $predicted_ORFs_bed;
my $predicted_ORFs_fasta;
my $plastid_image;

my $help = 0;

# Get command line arguments
GetOptions(
	'g=s'=>\$genome,
	'bam=s'=>\$bam_file,
	'db=s'=>\$blastdb,
	'sdir=s'=>\$script_dir,
	'wdir=s'=>\$workdir,
	'gtf=s'=>\$gtf,
	'en=s'=>\$experiment,
	'p=i'=>\$occupancy,
	'mn=i'=>\$min_read_len,
	'mx=i'=>\$max_read_len,
	'mo=i'=>\$MINORF,
	'mr=i'=>\$MINREAD,
	'ost=i'=>\$OFFSET_START,
	'osp=i'=>\$OFFSET_STOP,
	'osd=i'=>\$OFFSET_SD,
	'seed=s'=>\$SEED,
	'sd=s'=>\$USESD,
	'id=f'=>\$identity,
	'ev=f'=>\$evalue,
	'pg=i'=>\$pgm,
	'cdn=s'=>\$start_codons,
	'ncdn=s'=>\$start_codon_nset,
	'pcdn=s'=>\$start_codon_pset,
	'bgS=s'=>\$bedgraphS,
	'bgAS=s'=>\$bedgraphAS,
	'orf=s'=>\$predicted_ORFs,
	'bed=s'=>\$predicted_ORFs_bed,
	'fa=s'=>\$predicted_ORFs_fasta,
	'ps=s'=>\$plastid_image,
	'gcode=i'=>\$genetic_code,
  'by=s'=>\$seedBYpass,
  'score=f'=>\$score,
	'out=s'=>\$output_folder,
	'threads=i'=>\$threads,
	'h|help|?'=>\$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 1) if $help;

#----------------------------------------------------
#	Evaluate input variables
#----------------------------------------------------

# check mandatory variables
my %params = (
	g=>$genome,
	bam=>$bam_file,
	db=>$blastdb,
	sdir=>$script_dir
);


my @invalid = grep uninitialized_param($params{$_}), keys %params;
die "Not properly initialized: @invalid\n" if @invalid;


# check if mandatory variables are truly files
unless (-e $genome) {
    print "'$genome' is not a file. Please ensure the file exist\n";
    exit(1);
}

unless (-e $bam_file) {
    print "'$bam_file' is not a file. Please ensure the file exist\n";
    exit(1);
}

unless (-e $blastdb) {
    print "'$blastdb' is not a file. Please ensure the file exist\n";
    exit(1);
}

if ($gtf) {
    unless (-e $gtf) {
        print "'$gtf' is not a file. Please ensure the file exist\n";
        exit(1);
    }
}

unless ($occupancy == 1 or $occupancy == 3 or $occupancy == 5) {
	print "Option -p not properly initialize. -pg requires either 1, 3 or 5.\n";
    exit(1);
}

unless ($pgm == 1 or $pgm == 2) {
	print "Option -pg not properly initialize. -pg requires either 1 or 2.\n";
    exit(1);
}

if (uc($SEED) =~ m/^actgACTG/) {
	print "Option -seed not properly initialize. -seed string should contain only charatcers A,C,T, G.\n";
    exit(1);
}

unless (uc($USESD) eq 'Y' or uc($USESD) eq 'N') {
	print "Option -sd not properly initialize. -sd requires either Y or N.\n";
    exit(1);
}

unless ( $genetic_code >= 1 and $genetic_code <= 25) {
	print "Option -gcode not properly initialize. -gcode can take integer values between 1 and 25 (inclusive).\n";
    exit(1);
}

unless (uc($seedBYpass) eq 'Y' or uc($seedBYpass) eq 'N') {
	print "Option -by not properly initialize. -by requires either Y or N.\n";
    exit(1);
}


# check if script directoryis properly initialized and contains all scripts
if ($script_dir) {
	unless (-e $script_dir."/positive_set.pl") {
		print "Script 'positive_set.pl' not found in script directory.\nEnsure the directory  '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."/Ribo_seq_occupancy.py") {
		print "Script 'Ribo_seq_occupancy.py' not found in script directory.\nEnsure the directory  '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."/generate_all_ORFs.pl") {
		print "Script 'generate_all_ORFs.pl' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."/profiles.pl") {
		print "Script 'profiles.pl' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."/plot_profile.R") {
		print "Script 'plot_profile.R' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."/Random_Forest.R") {
		print "Script 'Random_Forest.R' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	unless (-e $script_dir."/post_processing.pl") {
		print "Script 'post_processing.pl' not found in script directory.\nEnsure the directory '$script_dir' contains all required file (see readme).\n";
		exit(1);
	}

	# Check prerequisits
	if ($pgm == 1) {
		check_if_pgm_exist('prodigal');
	} elsif ($pgm == 2) {
		check_if_pgm_exist('glimmer3');
	}

} else {
	print "Script directory not properly initialized. Please ensure the script directory is properly defined\n";
	exit(1);
}

# Generate bam index files
# Check for index file
if ((!-f $bam_file.".bai") || (-z $bam_file.".bai")) {
		print "No BAM index found. Computing... \n";
		my $command_index = "samtools index ".$bam_file;
		system($command_index) == 0
				or die ("Error running samtools: $! \n");
} else {
		print "BAM index file located. Advancing...\n";
}

# Prepare input variables
$experiment = ($experiment) ? $experiment."_": "";

$seedBYpass = uc($seedBYpass);
$USESD = uc($USESD);
$start_codons = uc($start_codons);	# convert to uppercase
$start_codon_nset = uc($start_codon_nset); # convert to uppercase

unless($start_codon_pset) {$start_codon_pset =$start_codons}
$start_codon_pset = uc($start_codon_pset);


# Check start codons input
my $start_cdn = {};
my @scodons = split /,/, $start_codons;
foreach my $codon(@scodons) {
    if (length($codon) != 3 or !(exists $translationHash{$codon})) {
        print "Codon '$codon' not a valid start codon\n";
	    print "Start codons must be 3 nucleotides long and contain either of A,C,G or T [example: -cdn ATG,GTG,TTG]\n";
        exit(1);
    }

	$start_cdn->{$codon} = 1;
}


# Check start codons for positive set
my $positive_codons = {};
my @pcodons = split /,/, $start_codon_pset;
foreach my $codon(@pcodons) {
    if (length($codon) != 3 or !(exists $translationHash{$codon})) {
        print "Codon '$codon' not a valid start codon for the positive set\n";
	    print "Positive set codon must be one or more of the standard Bacterial, Archaeal and Plant Plastid Code (transl_table=11 ATG,GTG,TTG).\n";
        exit(1);
    }
	$positive_codons->{$codon} = 1;
}


# Check start codons for negative set
my $negative_codons = {};
my @ncodons = split /,/, $start_codon_nset;
foreach my $codon(@ncodons) {
    if (length($codon) != 3 or !(exists $translationHash{$codon})) {
        print "Codon '$codon' not a valid start codon for the negative set\n";
	    print "Codons must be 3 nucleotides long and contain either of A,C,G or T [example: -ncdn ATG,GTG,TTG]\n";
        exit(1);
    }
	$negative_codons->{$codon} = 1;
}

# Check if working directory exist
my ($work_dir, $tmp_dir) = check_working_dir($workdir);

# append work_dir to output files
unless($bedgraphS) {$bedgraphS = $work_dir."/".$experiment."Ribo-seq_Sense_".$occupancy.".bedgraph";}
unless($bedgraphAS) {$bedgraphAS = $work_dir."/".$experiment."Ribo-seq_AntiSense_".$occupancy.".bedgraph";}
unless($predicted_ORFs) {$predicted_ORFs= $work_dir."/".$experiment."Predicted_ORFs.txt";}
unless($predicted_ORFs_bed) {$predicted_ORFs_bed = $work_dir."/".$experiment."Predicted_ORFs.bed";}
unless($predicted_ORFs_fasta) {$predicted_ORFs_fasta = $work_dir."/".$experiment."Predicted_ORFs.fasta";}
unless($plastid_image) {$plastid_image = $work_dir."/".$experiment."p_site_offset.png";}


# generate positive set
my $positive_set = $work_dir."/tmp/positive_set.txt";
print "Generate positive set...\n";
my $positive_set_gtf = $work_dir.'/tmp/positive.gtf';
my $cmd_positive = "perl ".$script_dir."/positive_set.pl $genome $blastdb $positive_set $min_read_len $max_read_len $MINORF $identity $evalue $start_codon_pset $pgm $work_dir $script_dir $threads $seedBYpass $genetic_code";
print "$cmd_positive\n\n";
system($cmd_positive);
print "Done.\n\n";


# section to implement plastid
my $psite_offset_file = $work_dir.$experiment."p_offsets.txt";
if ($occupancy == 1) {
	# check if plastid is installed
	my $search_psite = `which psite 2>&1`;
	chomp($search_psite);
	if ($search_psite =~ /^which: no/) {
		print "Could not locate ' psite '. Please ensure the plastid python package is installed.\n";
		exit(1);
	}

	my $search_metagene = `which metagene 2>&1`;
	chomp($search_metagene);
	if ($search_metagene =~ /^which: no/) {
		print "Could not locate ' metagene '. Please ensure the plastid python package is installed.\n";
		exit(1);
	}

	# find p-site offsets
    $psite_offset_file = generate_p_site($positive_set_gtf,$bam_file,$min_read_len,$max_read_len,$bam_file);
}


# Generate occupancy file
print "Generating ribosome occupancy file..\n";
my $occupancyFile = $work_dir."/".$experiment."Ribo-seq_".$occupancy."_occupancy.txt";
my $bedgraphS_prefix = $experiment."Ribo-seq_Sense_".$occupancy;
my $bedgraphAS_prefix = $experiment."Ribo-seq_AntiSense_".$occupancy;
my $cmd_occupancy = "python ".$script_dir."/Ribo_seq_occupancy.py $bam_file $occupancy $min_read_len $max_read_len $bedgraphS $bedgraphAS $occupancyFile $bedgraphS_prefix $bedgraphAS_prefix $psite_offset_file";
print "$cmd_occupancy\n\n";
system($cmd_occupancy);
unless (-e $occupancyFile) {
    print "'$occupancyFile' file does not exist.\n";
    exit(1);
}
print "Done.\n\n";


# Generate all possible ORFs
print "Generate all possible ORFs...\n";
my $codons = $start_codons.",".$start_codon_pset.",".$start_codon_nset;	# combine the start codons
my $ORF_file = $work_dir."/tmp/all_valid_ORFs.txt";
my $cmd_orf_gen = "perl ".$script_dir."/generate_all_ORFs.pl $genome $occupancyFile $positive_set $ORF_file $MINORF $OFFSET_START $OFFSET_STOP $OFFSET_SD $SEED $codons $work_dir $script_dir $threads";
print "$cmd_orf_gen\n\n";
print "\n";
system($cmd_orf_gen);
print "Done.\n\n";


# Generate metagenic profile
print "Meta-genic plots\n";
my $cmd_meta = "perl ".$script_dir."/profiles.pl $positive_set $occupancyFile $MINREAD $work_dir $script_dir";
print "$cmd_meta\n\n";
print "\n";
system($cmd_meta);
print "Done.\n\n";


# ORF prediction
print "Performing ORF prediction analysis..\n";
my $RF_prediction = $work_dir."/tmp/RF_predicted_ORFs.txt";
my $threshold = $work_dir."/tmp/threshold.txt";

my $RF_command = "Rscript ".$script_dir."/Random_Forest.R $ORF_file $positive_set $work_dir $start_codons $start_codon_nset $USESD $score";
print "$RF_command\n\n";
print "\n";
system($RF_command);
unless (-e $RF_prediction) {
    print "'$RF_prediction' file does not exist.\n";
    exit(1);
}
print "Done.\n\n";


# Post processing
print "Cleaning up RF predictions..\n";
my $output_prefix = $work_dir."/".$experiment."Predicted_ORFs";
my $processing_cmd = "";
if ($gtf) {
    $processing_cmd = "perl ".$script_dir."/post_processing.pl $RF_prediction $genome $occupancyFile $output_prefix $threshold $OFFSET_START $MINREAD $predicted_ORFs $predicted_ORFs_bed $predicted_ORFs_fasta $gtf";
} else {
    $processing_cmd = "perl ".$script_dir."/post_processing.pl $RF_prediction $genome $occupancyFile $output_prefix $threshold $OFFSET_START $MINREAD $predicted_ORFs $predicted_ORFs_bed $predicted_ORFs_fasta";
}
print "$processing_cmd\n\n";
system($processing_cmd);
print "Done.\n\n";


timer($startRun);	# Get total Run time


##################
##	SUBS
##################

## Generate metagene and p-site estimates
sub generate_p_site {
	  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
		my $logs_dir = $work_dir."/logs";
		system("mkdir -p $logs_dir");

    my $genes_gtf = $_[0];
    my $bam = $_[1];
    my $min_l = $_[2];
    my $max_l =$_[3];
		my $bam_file =$_[4];

    # temporary files
    my $run_name = $work_dir."/tmp/plastid";

		my $sort_tmp_file = $work_dir."/tmp/check_sort.txt";

		my $log_metagene = $work_dir."/logs/metagene$min$hour$mday$mon$year.log";
    #Build command
    my $command_meta = "metagene generate -q ".$run_name." --landmark cds_start --annotation_files ".$genes_gtf." 2> $log_metagene";
    print "$command_meta\n";
    system($command_meta) == 0
        or die ("Error running metagene: $! \n");

    my $log_psite = $work_dir."/logs/psite$min$hour$mday$mon$year.log";
    #Build command
    my $psitefile = $run_name."_rois.txt";
    my $command_psite = "psite ".$run_name."_rois.txt ".$run_name." --min_length ".$min_l." --max_length ".$max_l." --require_upstream --count_files ".$bam_file." 2> $log_psite";
    print "$command_psite\n";
    my $psite_off_output = $work_dir."/".$experiment."p_site_offsets.txt";
    my $offsetdefaults = "length\tp_offset\n22\t19\n23\t20\n24\t21\n25\t22\n26\t23\n27\t24\n28\t13\n29\t13\n30\t13\n31\t13\n32\t13\n33\t13\n34\t13\n35\t13\n36\t13\n37\t13\n38\t21\n39\t22\n40\t23\ndefault\t13\n";
    my $offsetfile = $run_name. "_p_offsets.txt";
    unless(system($command_psite)  == 0){
        open( my $fh, '>', $offsetfile) or die "Could not open file $offsetfile $!";
        print $fh $offsetdefaults;
        close $fh;
        system("touch " . $run_name . "_p_offsets.png;");
    }
        #or die ("Error running psite: $! \n");
    print "touch " . $run_name. "_p_offsets.txt; touch " . $run_name . "_p_offsets.png;";
    system("cp " . $run_name . "_p_offsets.txt $psite_off_output") or print ("Error copying psite_off_output: " . $run_name . "_p_offsets.txt $psite_off_output $!");
    system("cp " . $run_name . "_p_offsets.png $plastid_image") or print ("Error copying p_offset.png: " . $run_name . "_p_offsets.png $plastid_image $!");
    return $psite_off_output;
}


sub check_if_pgm_exist {

	my $pgm = $_[0];

	# check for prodigal or glimmer
	my $search = `which $pgm 2>&1`;
    print "$search\n";
	chomp($search);
	if ($search =~ /^which: no/) {
		if ($pgm eq "prodigal") {
			unless (-x $script_dir."/bin/prodigal") {	# if prodigal not install
				print "Could not locate ' $pgm '. Please ensure the program is installed or present in the script directory and it is executable: $! \n";
				exit(1);
			}
		} elsif ($pgm eq "glimmer3") {
            print "$script_dir"."/bin/glimmer/glimmer3\n";
			unless (-e $script_dir."/bin/glimmer/glimmer3" and -e $script_dir."/bin/glimmer/build-icm") {
				print "Could not locate ' $pgm '. Please ensure the program is installed or present in the script directory and it is executable: $! \n";
				exit(1);
			}
		}
	}
}


sub check_working_dir {

	my $work_dir = $_[0];
	my $tmp_dir;

	if ($work_dir) {
		$tmp_dir = $work_dir."/tmp";
		if (!-d $work_dir) {
			system("mkdir -p $work_dir" or die "Couldn't create '$work_dir.\n");
			system("mkdir -p $tmp_dir" or die "Couldn't create '$tmp_dir'.\n");
		} else {
			system("rm -rf $work_dir" or die "Can delete '$work_dir': $!\n");
			system("mkdir -p $work_dir" or die "Couldn't create '$tmp_dir'.\n");
			system("mkdir -p $tmp_dir" or die "Couldn't create '$tmp_dir'.\n");
		}
	} else {
		$work_dir = getcwd();
		$work_dir = $work_dir."/".$experiment.$output_folder;
		$tmp_dir = $work_dir."/tmp";

		if (!-d $work_dir) {
			system("mkdir -p $work_dir");
			system("mkdir -p $tmp_dir");
		}
	}

	$work_dir = $work_dir."/";
	$tmp_dir = $tmp_dir."/";
	return $work_dir, $tmp_dir;
}


sub uninitialized_param {
	my ($v) = @_;
	not ( defined($v) and length $v );
}

sub timer {
	my $startRun = shift;
	my $endRun 	= time();
	my $runTime = $endRun - $startRun;
	printf("\nTotal running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime  % 3600) / 60), int($runTime % 60));
}

__END__

=head1 NAME

reparation - Ribosome Profiling Assisted (Re-)Annotation of Bacterial genomes.

=head1 SYNOPSIS

reparation.pl [options] -bam riboseq_alignment_files_bam_format -g bacteria_genome_fasta_file -sdir scripts_directory -db curated_protein_db(fasta)

=head1 OPTIONS

=over

=head2 Mandatory

=item B<-bam>

Ribosome alignment file (bam)

=item B<-g>

Genome fasta file. This should be the same genome fasta file used in the alignment of the Ribo-seq reads.

=item B<-sdir>

The "scripts" directory (avialable within the REPARATION directory), defaults to directory of reparation.pl script

=item B<-db>

fasta database of curated bacteria protein sequences

=head2 Optional

=item B<-gtf>

GTF genome annotation file

=item B<-wdir>

working directory (defaults to current directory)

=item B<-en>

Experiment name

=item B<-p>

Ribosome profiling read p site assignment strategy, 1 = plastid P-site estimation ((default), 3 = 3 prime of read, 5 prime of the read

=item B<-mn>

All ribosome profiling reads shorter than these values are eliminated from the ananlysis (default = 22)

=item B<-mx>

All ribosome profiling reads longerer than these values are eliminated from the ananlysis (default = 40)

=item B<-mo>

Minimum length of open reading frame considered for prediction (default = 30 value in nucleotides)

=item B<-mr>

Open reading frames with less than these number of ribosome profiling reads are eliminated from analysis (default = 3)

=item B<-ost>

Start region length in nucleotides (default = 45nts). This value is used to calculate features specific to the start region.

=item B<-osp>

Stop region length in nucleotides (default = 21nts). This value is used to calculate features specific to the stop region.

=item B<-osd>

Distance of Shine dalgarno sequence from start codon (defualt = 5nts).

=item B<-seed>

Shine dalgarno sequence (default = GGAGG). The shine dalgarno sequence used for Ribosome binind site energy calculation.

=item B<-sd>

Use ribosome binding site (RBS) energy in the open reading frame prediction (Y = use RBS energy (default), N = donot use RBS engergy)

=item B<-id>
Online-Ticket
UNI STUTT
Minimum identify score for BLAST protein sequence search (default = 0.75)

=item B<-ev>

maximum e-vlaue for BLAST protein sequence search (default = 1e-5)

=item B<-pg>

program for initial positive set generation (1 = prodigal (default), 2 = glimmer)

=item B<-cdn>

Comma separted list of start codons (default = ATG,GTG,TTG)

=item B<-ncdn>

Start codon for negative set (default = CTG)

=item B<-pcdn>

Start codon for positive set (default = ATG,GTG,TTG). Should be a subset of the standard genetic code for bacteria

=item B<-gcode>

Genetic code to use for initail positive set generation. Valid when -pg is 1. (default = 11, takes value between 1-25)

=item B<-by>

Flag to determine if prodigal should bypass Shine-Dalgarno trainer and force a full motif scan (default = N). (Y = yes, N = no) Valid only for -pg 1

=item B<-score>

Random forest classification probability score threshold to define as ORF are protein coding, the minimum  (defualt is 0.5)

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut
