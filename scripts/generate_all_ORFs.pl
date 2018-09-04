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

use Bio::SeqIO;
use Getopt::Long;
use POSIX;

# Get command line arguments
my ($genome,$occupancyFile,$positive_set,$ORF_file,$MINORF,$OFFSET_START,$OFFSET_STOP,$OFFSET,$SEED,$codons,$work_dir,$script_dir, $threads) = @ARGV;


#----------------------------------------------------
#			EXECUTION
#----------------------------------------------------

my %translationHash = 
	(GCA => "A", GCG => "A", GCT => "A", GCC => "A", GCN => "A",
     TGC => "C", TGT => "C",
     GAT => "D", GAC => "D",
     GAA => "E", GAG => "E",
     TTT => "F", TTC => "F",
     GGA => "G", GGG => "G", GGC => "G", GGT => "G", GGN => "G",
     CAT => "H", CAC => "H",
     ATA => "I", ATT => "I", ATC => "I",
     AAA => "K", AAG => "K",
     CTA => "L", CTG => "L", CTT => "L", CTC => "L", CTN => "L", TTA => "L", TTG => "L",
     ATG => "M",
     AAT => "N", AAC => "N",
     CCA => "P", CCT => "P", CCG => "P", CCC => "P", CCN => "P",
     CAA => "Q", CAG => "Q",
     CGA => "R", CGG => "R", CGC => "R", CGT => "R", CGN => "R",
     AGA => "R", AGG => "R",
     TCA => "S", TCG => "S", TCC => "S", TCT => "S", TCN => "S",
     AGC => "S", AGT => "S",
     ACA => "T", ACG => "T", ACC => "T", ACT => "T", ACN => "T",
     GTA => "V", GTG => "V", GTC => "V", GTT => "V", GTN => "V",
     TGG => "W",
     TAT => "Y", TAC => "Y",
     TAG => "*", TAA => "*", TGA => "*");

# Stop codons
my $stop_codons = {TAG => "*", TAA => "*", TGA => "*"};

# Process start codons
my $start_codons = {};
my @codons = split /,/, $codons;
foreach my $codon(@codons) {
	$start_codons->{$codon} = 1;
}

# check if input files were properly generated
unless (-e $positive_set) {
    print "'$positive_set' file does not exist.\n";
    exit(1);
}

unless (-e $occupancyFile) {
    print "'$occupancyFile' file does not exist.\n";
    exit(1);
}

# read genome fasta file
my $genomes = read_fasta($genome);

# read RPF positional information from file
my ($reads_table,$mapped_total) = get_read_table($occupancyFile);


# generate positive set 
my $ORF_positive = get_positive_set($positive_set);

# Calculate RBS energy
my $WINDOW 		= 20 - $OFFSET;		# the window size to look RBS for
my $SEEDLENGTH 	= length($SEED);
my $THRESHOLD 	= 3*($SEEDLENGTH - 2);
my ($SD_markov_score, $Max_loc_score, $markov_score_threshold, $perfect_markov_score, $initial_threshold) = Train_SD_Markov_Model($ORF_positive);

print "Markov Threshold score $markov_score_threshold\n";
print "Prefect Markov Score $perfect_markov_score\n";

# generate all ORFs
my $max_length = 1;
my $ORFs = find_all_ORFs($genomes);

write_ORFs_to_file($ORFs, $ORF_file);

print "ORF Ribo-seq info written to file.\n";


##################
##	SUBS
##################


sub write_ORFs_to_file {

	my $ORFs 	 = $_[0];
	my $filename = $_[1];

	my $count_orf = 0;
	my $count_nterm = 0;

	open(F, ">$filename") or die ("Can't creat file $filename: $!\n");
	print F "orf_id\tgene\tstrand\tlength\tstart_codon\tcount\trpkm\tcoverage\tstart_coverage\tstart_rpkm\tstop_rpkm\taccumulation_proportion\tSD_score\tSD_pos\n";
	foreach my $ORF (sort keys %$ORFs) {
		my $length 		= $ORFs->{$ORF}->{len};
		my $rpkm 		= $ORFs->{$ORF}->{rpkm};
		my $region 		= $ORFs->{$ORF}->{region};
		my $strand 		= $ORFs->{$ORF}->{strand};
		my $start_codon	= $ORFs->{$ORF}->{start_codon};
		my $start 		= $ORFs->{$ORF}->{start};
		my $stop 		= $ORFs->{$ORF}->{stop};
		my $coverage 	= $ORFs->{$ORF}->{coverage};
		my $count 		= $ORFs->{$ORF}->{count};
		my $start_rpkm	= $ORFs->{$ORF}->{start_rpkm};
		my $stop_rpkm	= $ORFs->{$ORF}->{stop_rpkm};
		my $start_cov	= $ORFs->{$ORF}->{start_coverage};
		my $start_count	= $ORFs->{$ORF}->{start_count};
		my $gene 		= $ORFs->{$ORF}->{gene};

		# calculate RBS energy
		my ($SD_seq,$SD_score,$SD_pos) = calculate_RBS($ORF,$strand);
		my ($accumu_prop,$average_start, $average_rest) = Accumulation_proportion($start,$stop,$strand,$region,$count);

		print F "$ORF\t$gene\t$strand\t$length\t$start_codon\t$count\t$rpkm\t$coverage\t$start_cov\t$start_rpkm\t$stop_rpkm\t$accumu_prop\t$SD_score\t$SD_pos\n";

		$count_orf++;
	}
	close F;

	print "Total number of ORFs generated $count_orf.\n";

}


sub Accumulation_proportion {

	my $start = $_[0];
	my $stop = $_[1];
	my $strand = $_[2];
	my $region = $_[3];
	my $count = $_[4];

	my $length = $stop - $start + 1;
	my $tmp_OFFSET_START = $OFFSET_START;
	if ($length < $OFFSET_START) {$tmp_OFFSET_START = ceil(0.7*$length)}
	my $start_start = ($strand eq '+') ? $start: $stop - $tmp_OFFSET_START;
	my $start_stop = ($strand eq '+') ? $start + $tmp_OFFSET_START: $stop;
	my $rest_start = ($strand eq '+') ? $start + $tmp_OFFSET_START + 1: $start;
	my $rest_stop = ($strand eq '+') ? $stop : $stop - $tmp_OFFSET_START - 1;

    my $orf_average_reads = $count/$length;

	my $count_rest = 0;
	my $count_start = 0;
	for ( my $pos = $start; $pos <= $stop; $pos++) {
		if (exists $reads_table->{$region}->{$strand}->{$pos}) {
			if ($start_start <= $pos and $pos <= $start_stop) {
                my $reads = $reads_table->{$region}->{$strand}->{$pos};
				$count_start += $reads;
			}

			if ($rest_start <= $pos and $pos <= $rest_stop) {
                my $reads = $reads_table->{$region}->{$strand}->{$pos};
                #if ($reads >= 0.7*$count) {$reads = $orf_average_reads}
				$count_rest += $reads;
			}
		}
	}

	my $average_start = $count_start/(abs($start_stop - $start_start) + 1);
	my $average_rest = $count_rest/(abs($rest_stop - $rest_start) + 1);

	if ($average_rest <= 1 or $average_start<=1) {
		return 0, $average_start, $average_rest;
	} else {
		return $average_start/$average_rest, $average_start, $average_rest;
	}

}


sub find_all_ORFs {
	
	my $genome 	= $_[0];
    my $ORFs = {};

	foreach my $region (keys %$genome) {
		my $directStrand = $genome->{$region};
		my $reverseComplement = revdnacomp($directStrand);

		for (my $i = 0; $i < 2; $i = $i + 1) {
			my $sequenceEntry = "";
			my $strand = "";
			if ($i == 0) {		# forward strand.
				$sequenceEntry = $directStrand;
				$strand = "+";
			} else {			# reverse strand
				$sequenceEntry = $reverseComplement;
				$strand = "-";
			}

			my @starts;
			foreach my $start_cdn (keys %$start_codons) {
				my $offset = 0;
				my $start = index($sequenceEntry, $start_cdn, $offset);
				push @starts, $start;
				while ($start != -1) {
					$offset = $start + 1;
					$start = index($sequenceEntry, $start_cdn, $offset);
					push @starts, $start;
				}
			}

			# check if there is a corresponding stop to all start positions
			foreach my $start (@starts) {
				my $stop;
				my $aa_seq = "";		# Amino acide sequence
				my $tr_seq = "";		# DNA sequence
				my $strt_cdn = uc(substr($sequenceEntry, $start, 3));

				for (my $i = $start; $i <= (length($sequenceEntry) - 3); $i = $i + 3) {
					my $codon = uc(substr($sequenceEntry, $i, 3));
					last if (length($codon) < 3);
					if (!(exists $translationHash{$codon})) {
                        $aa_seq = "";
                        $tr_seq = ""; 
                        last;
                    }

					my $aa = $translationHash{$codon};
					if ($aa eq "*") {
						$stop = $i;
						last;
					}
					$aa_seq = $aa_seq.$aa;
					$tr_seq = $tr_seq.$codon;
				}

				if ($stop) {	# keep tract only on start with corresponding stop sequences
					my $length = abs($stop - $start);
					my $new_start = $start;
					my $new_stop = $stop;
					if ($length >= $MINORF) {
						if ($strand eq '-') {
							my $start_tmp = length($directStrand) - $new_stop + 1;
							$new_stop  = length($directStrand) - $new_start;
							$new_start = $start_tmp;
						} else {
							$new_start += 1;
						}

						my ($count,$coverage,$rpkm) = RIBO_info($region, $strand, $new_start, $new_stop,1);

						my $tmp_OFFSET_START = $OFFSET_START;
						my $tmp_OFFSET_STOP = $OFFSET_STOP;
						if ($length < $OFFSET_START + $OFFSET_STOP) {
							$tmp_OFFSET_START = ceil(0.7*$length);
							$tmp_OFFSET_STOP = ceil(0.25*$length);
						}
		
						# get start region interval adjust per strand
						my $start_start = $new_start;
						my $start_stop =  $new_stop;
						if ($strand eq '-') {$start_start = $new_stop - $tmp_OFFSET_START}
						if ($strand eq '+') {$start_stop = $new_start + $tmp_OFFSET_START}
						my ($start_count,$start_coverage, $start_rpkm) = RIBO_info($region, $strand, $start_start, $start_stop,$count);

						# get stop region interval adjust per strand
						my $stop_start = $new_start;
						my $stop_stop = $new_stop;
						if ($strand eq '+') {$stop_start = $new_stop - $tmp_OFFSET_STOP} 
						if ($strand eq '-') {$stop_stop =  $new_start + $tmp_OFFSET_STOP}
						my ($stop_count,$stop_coverage, $stop_rpkm) = RIBO_info($region, $strand, $stop_start, $stop_stop,$count);

						my $id = $region.":".$new_start."-".$new_stop;	# ORF location id
						$ORFs->{$id}->{region} 		= $region;
						$ORFs->{$id}->{start} 		= $new_start;
						$ORFs->{$id}->{stop} 		= $new_stop;
						$ORFs->{$id}->{len} 		= $length;
						$ORFs->{$id}->{strand} 		= $strand;
						$ORFs->{$id}->{rpkm} 		= $rpkm;
						$ORFs->{$id}->{start_rpkm} 	= $start_rpkm;
						$ORFs->{$id}->{stop_rpkm} 	= $stop_rpkm;
						$ORFs->{$id}->{count} 		= $count;
						$ORFs->{$id}->{start_count} = $start_count;
						$ORFs->{$id}->{aa_seq} 		= $aa_seq;
						$ORFs->{$id}->{coverage} 	= $coverage;
						$ORFs->{$id}->{start_coverage} = $start_coverage;
						$ORFs->{$id}->{start_codon} = $strt_cdn;
						$ORFs->{$id}->{gene} 		= ($strand eq "+") ? $region.":FAM+".$new_stop: $region.":FAM-".$new_start;
						if ($max_length < $length) {$max_length = $length;}
					}
				}
			}
		}
	}

	return $ORFs;
}



sub RIBO_info {
	# Function to read RIBO-seq positional data into memory
	my $chr 	  = $_[0];
	my $strand 	  = $_[1];
	my $start 	  = $_[2];
	my $stop 	  = $_[3];
	my $orf_count = $_[4];

	my $count = 0;
	my $coverage = 0;
	my $rpkm = 0;

	if ($strand eq '+') {
		$start -= 3;
	} elsif ($strand eq '-') {
		$stop += 3;
	}

	if ($orf_count == 0) {
		return $count,$coverage,$rpkm;
	} else {
		for ( my $pos = $start; $pos <= $stop; $pos++) {
			if (exists $reads_table->{$chr}->{$strand}->{$pos}) {
				$count += $reads_table->{$chr}->{$strand}->{$pos};
				$coverage++;
			}
		}
		my $length = $stop - $start + 1;
		$coverage 	= $coverage/$length;
		$rpkm = (($count/$orf_count)*1000000000)/($length*$mapped_total);

		return $count,$coverage,$rpkm;
	}
}

sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}


sub get_read_table {
	my $file = $_[0];
	my $read = {};		# hash to store all ribo-seq inofrmation from SQLite DB
	my $mapped_total;

	open(F, "$file") or die "Cannot open file $file\n";
	while (<F>) {
		chomp $_;
		if ( $_ =~ /^mappable/) {
			$mapped_total = (split ':', $_)[1];
		} else {
			my @line = split '\t', $_;
			my $region 	= $line[0];
			my $start 	= $line[1];
			my $strand 	= $line[2];
			$read->{$region}->{$strand}->{$start} = $line[3];
		}
	}
	
	return $read, $mapped_total;
}


sub read_fasta {
	my $file = $_[0];
	my $cdna = {};

	my $in  = Bio::SeqIO->new(-file => $file, -format => "fasta");
	while(my $seqs = $in->next_seq) {
		my $id  = $seqs->display_id;	
		my $seq = $seqs->seq;
		my $desc = $seqs->desc;

		$cdna->{$id} = uc($seq);
	}
	return $cdna;
}


#-------------- Positive Set -------------#


sub get_positive_set {

	my $positive_set = $_[0];
	my $ORF_positive = {};

	open (my $F, $positive_set) or die  "Error reading file: $positive_set\n";
	while (<$F>) {
		next if (/^orf_id/);
		chomp $_;
		my @line = split '\t', $_;
		my $orf_id = $line[0];
		my $strand = $line[1];

		$ORF_positive->{$orf_id} = $strand;
	}
	close $F;
	return $ORF_positive;
}



#-------------- Calculate RBS energy -------------#

sub calculate_RBS {

	my $ORF 	= $_[0];
	my $strand	= $_[1];

	my ($start)	= $ORF =~ /:(\d+)-/;
	my ($stop)	= $ORF =~ /-(\d+)$/;
	my ($region)= $ORF =~ /^(.*):/;

	my $SD_region =  ($strand eq '+') ? get_SD_region($start, $strand, $region): get_SD_region($stop, $strand, $region);

	my ($SD_seq, $SD_score, $SD_pos);

	if (defined $SD_region and length($SD_region) >= $SEEDLENGTH) {
		($SD_seq, $SD_score, $SD_pos) = calculate_SD_score($SD_region);
		$SD_pos = ($strand eq '+') ? $start - $SD_pos - $OFFSET: $stop + $OFFSET + $SD_pos;
	} else {
		($SD_seq, $SD_score, $SD_pos) = ("-",-100000000,'NA');
	}

	return $SD_seq,$SD_score,$SD_pos;
}


sub calculate_SD_score {

	my $SD_region = $_[0];

	my $kmers = kmers($SD_region);	# get all kmers within the SD region
	my $max_SD_score;
	my $max_SD_pos;
	my $max_SD_seq;

	foreach my $p (keys %$kmers) {
		my $kmer_seq = $kmers->{$p};
		my $kmer_SD_score = SD_kmer_score($kmer_seq, $p);

		unless (defined $max_SD_seq) {		# if first score
			$max_SD_pos	 = $p;
			$max_SD_score= $kmer_SD_score;
			$max_SD_seq	 = $kmer_seq;
		} else {
			if ($max_SD_score < $kmer_SD_score) {
				$max_SD_pos	 = $p;
				$max_SD_score= $kmer_SD_score;
				$max_SD_seq	 = $kmer_seq;
			}
		}
	}

	return $max_SD_seq, $max_SD_score, $max_SD_pos;
}


sub SD_kmer_score {
	my $seq = $_[0];
	my $pos = $_[1];

	my @nucl = split '', $seq;
	my $scorei;
	for(my $i=0; $i < scalar(@nucl); $i++) {
		if ($scorei) {
			if ($SD_markov_score->{$i}->{$nucl[$i]}) {
                $scorei = $scorei*$SD_markov_score->{$i}->{$nucl[$i]};
            }
		} else {
			$scorei = $SD_markov_score->{$i}->{$nucl[$i]};
        }
	}

	return($Max_loc_score->{$pos} + $scorei);
}


sub Train_SD_Markov_Model {

	my $positive_ORFs = $_[0];

	my @sequences = get_training_sequences($positive_ORFs);

	my $Max_loc_score 	= {};
	my $SD_markov_score = {};
	my @training_seq;

	foreach my $seq (@sequences) {
		my $MAX_score =  -100000000;	# default maximum similarity score
		my $MAX_loc;					# position from start with maximum similarity score
		my $MAX_seq;					# Shine dalgano sequence with max similarity score

		my $kmers = kmers($seq);		# get all kmers within the SD region
		foreach my $p (keys %$kmers) {

			my $kmer_seq = $kmers->{$p};
			my $sim_score = simple_similarity($kmer_seq);

			if ($MAX_score < $sim_score) {
				$MAX_loc 	 = $p;
				$MAX_score 	 = $sim_score;
				$MAX_seq	 = $kmer_seq;
			}
		}

		# Collect all sequences scoring above threshold
		if ($MAX_score >= $THRESHOLD) {
			push @training_seq, $MAX_seq;
			$Max_loc_score->{$MAX_loc}++;	# get frequencies of each location in the training set
		}
	}

	# Calculate frequency for each nucleotide
	foreach (@training_seq) {
		my @SD_seq = split '',$_;
		for(my $i = 0; $i < scalar(@SD_seq); $i++) {
			$SD_markov_score->{$i}->{$SD_seq[$i]}++;
		}
	}

	my $N = scalar(@training_seq);

	#---- LOGARITHMIC SCORING ---------------------
	for( my $i=0; $i < $SEEDLENGTH; $i++) {
		$SD_markov_score->{$i}->{'A'} = log(($SD_markov_score->{$i}->{'A'}+0.000001)/$N);
		$SD_markov_score->{$i}->{'C'} = log(($SD_markov_score->{$i}->{'C'}+0.000001)/$N);
		$SD_markov_score->{$i}->{'G'} = log(($SD_markov_score->{$i}->{'G'}+0.000001)/$N);
		$SD_markov_score->{$i}->{'T'} = log(($SD_markov_score->{$i}->{'T'}+0.000001)/$N);
	}

	# Score the locations within the DS_region
	for( my $i=0; $i < $WINDOW; $i++) {
		unless ($Max_loc_score->{$i}) {$Max_loc_score->{$i}=0}
		$Max_loc_score->{$i} = log(($Max_loc_score->{$i}+0.0001)/$N);
	}

    # determine minimum markov score for threshold
    my $min_score = 100000000;
    my $max_score = -100000000000;
    my $skip_val=log(0.000001/$N);
	foreach my $p (keys %$SD_markov_score) {
		foreach my $nucl (keys %{$SD_markov_score->{$p}}) {
			
			next if ($SD_markov_score->{$p}->{$nucl} <= $skip_val); # skip if default score i.e not observed in training set
			if ($SD_markov_score->{$p}->{$nucl} > $max_score) { # maximum score
				$max_score = $SD_markov_score->{$p}->{$nucl};
			}
			
			if ($SD_markov_score->{$p}->{$nucl} < $min_score) { # minimun score
				$min_score = $SD_markov_score->{$p}->{$nucl};
			}
		}
	}

    my $min_loc_posib=100000000000;	
    my $max_loc_posib=-100000000000;

    $skip_val=log(0.0001/$N);
	foreach my $p (keys %$Max_loc_score) {
		next if ($Max_loc_score->{$p} <= $skip_val); 
		if ($Max_loc_score->{$p} > $max_loc_posib) { # maximum score
			$max_loc_posib = $Max_loc_score->{$p};
		}
			
		if ($Max_loc_score->{$p} < $min_loc_posib) { # minimun score
			$min_loc_posib = $Max_loc_score->{$p};
		}		
	}

    my $markov_score_threshold  = ($min_score+$max_score)*0.5+$min_loc_posib;
    my $perfect_markov_score 	= $max_score+$max_loc_posib;
    my $initial_threshold 		= ($perfect_markov_score+$markov_score_threshold)*0.5;

	print "Number of sequences used in SD scoring training set $N\n";

	return $SD_markov_score, $Max_loc_score, $markov_score_threshold, $perfect_markov_score, $initial_threshold;
}


sub get_training_sequences {
	# collect all SD region from coordinates in the training set;
	my $positive_ORFs  = $_[0];

	my @training_seq;
	foreach my $ORF (keys %$positive_ORFs) {

		my $strand = $positive_ORFs->{$ORF};
		unless ($strand) {print "$ORF\n";next}

		my ($start)= $ORF =~ /:(\d+)-/;
		my ($stop)= $ORF =~ /-(\d+)$/;
		my ($region)= $ORF =~ /^(.*):/;

		if ($strand eq '+') {
			my $context = get_SD_region($start, $strand, $region);
			next unless ($context);
			next if (length($context) < $SEEDLENGTH);
			push @training_seq, $context;
		} elsif ($strand eq '-') {
			my $context = get_SD_region($stop, $strand, $region);
			next unless ($context);
			next if (length($context) < $SEEDLENGTH);
			push @training_seq, $context;
		}
	}

	return @training_seq;
}


sub simple_similarity {

	# function to calculate the simple similarty score of a sequence
	# relative to the Shine dalgano seed sequence
	# the score is weighted based on the number of hydrogen bonds

	my $seq = $_[0];

	my $score = 0;
	for (my $i = 0; $i < length($seq); $i++) {
		if (substr($SEED, $i,1) eq substr($seq, $i,1)) {
			if (substr($SEED, $i,1) eq "G" or substr($SEED, $i,1) eq "C") {
				$score += 3;
			} elsif (substr($SEED, $i,1) eq "A" or substr($SEED, $i,1) eq "T") {
				$score += 2;
			}
		} else {
			if (substr($SEED, $i,1) eq 'A' and substr($seq, $i,1) eq 'G') {
				$score += 1.5;
			} elsif (substr($SEED, $i,1) eq 'C' and substr($seq, $i,1) eq 'T') {
				#$score += 1.5;
			}
		}
	}

	return $score;
}


sub kmers {
	# get all kmers in the SD region
    my $seq = $_[0];

    my $kmers;
    for (my $i = length($seq); $i >= $SEEDLENGTH; $i--) {
		last if ($i - $SEEDLENGTH < 0);
		my $pos = $WINDOW - $i;				# position closest to ORF start
        $kmers->{$pos} = substr($seq, $i-$SEEDLENGTH, $SEEDLENGTH);
    }

    return $kmers;
}


sub get_SD_region {

	# Function to get the Shine Dalgano sequnce region
	my $start	= $_[0];
	my $strand 	= $_[1];
	my $region 	= $_[2];

	my $context;
	if ($strand eq "+") {
		my $seq = $genomes->{$region};
		my $start_offset = $start - $OFFSET - $WINDOW;

		my $WINDOW_tmp = $WINDOW;
		if ($start_offset < 0) {
			$WINDOW_tmp = $WINDOW_tmp + $start_offset;	# adjust the length for the overhang
			$start_offset = 0;
		}
		if ($WINDOW_tmp > 0) { $context = substr($seq, $start_offset, $WINDOW_tmp);}

	} elsif ($strand eq "-") {
		my $seq = $genomes->{$region};
		my $start_offset = $start + $OFFSET;
		my $WINDOW_tmp = $WINDOW;
		if ($start_offset > length($seq)) {
			$WINDOW_tmp = $WINDOW_tmp + (length($seq) - $start_offset);
			$start_offset = 0;
		}
		if ($WINDOW_tmp > 0) {$context = revdnacomp(substr($seq, $start_offset, $WINDOW_tmp));}
	}

	return $context;

}




