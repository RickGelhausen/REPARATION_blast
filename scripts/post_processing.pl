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
use POSIX;


#############################################################################
#
#############################################################################


my ($RF_prediction,$genome,$occupancyFile,$output_prefix,$threshold,$OFFSET_START,$MINCOUNT, $predicted_ORFs, $predicted_ORFs_bed, $predicted_ORFs_fasta, $gtf_file) = @ARGV;


my $WINDOW = 21;	# Window to calculate kurtosis
my $REGION = 30;	# Region to account for any mapping errors we assume the aprox. length of a ribosome.
my $SDLEN = 10;		# Approximate length of the Shine dalgarno sequence
my $OF_FOLD = 5;

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
     TAG => "", TAA => "", TGA => "");

# read genome and RPF occupancy files
my $genomes = read_fasta($genome);
my ($reads_table,$mapped_total) = get_read_table($occupancyFile);

# Get thresholds from file
open(F, $threshold) or die "Cannot open file $threshold\n";
my @thresholds = <F>;
my $MINCOV = $thresholds[0];
my $MINRPKM = $thresholds[1];

# Random Forest ORFS
my $ORFs_RF = read_file($RF_prediction);

my $ORFs_processed = ORFs_start_selection($ORFs_RF);

$ORFs_processed = internal_out_of_frame($ORFs_processed);
$ORFs_processed = two_way_overlap($ORFs_processed);
$ORFs_processed = one_way_overlap($ORFs_processed);

if ($gtf_file) {	# check if gtf is a file

	my $annotated = annotation($gtf_file);
	my $Classified_ORFs = ORF_classification($ORFs_processed, $annotated);

	create_bed_classified($Classified_ORFs, $output_prefix);
	output_ORFs_classified($Classified_ORFs, $output_prefix);
	predicted_ORFs_pep($ORFs_processed,$output_prefix);

	my $count_Annotated = 0;
	my $count_Novel = 0;
	my $count_Intergenic = 0;
	my $count_Internal = 0;
	my $count_Truncation = 0;
	my $count_ncRNA = 0;
	my $count_Extension = 0;
	my $count_Reverse = 0;
	my $count_3_codon = 0;
	my $count_out_of_frame = 0;

	foreach my $gene (keys %$Classified_ORFs) {
		foreach my $ORF (keys %{$Classified_ORFs->{$gene}}) {
			if ($Classified_ORFs->{$gene}->{$ORF}->{orf_type} eq 'Annotated') {
				$count_Annotated++;
			} elsif ($Classified_ORFs->{$gene}->{$ORF}->{orf_type} eq 'Truncation') {
				if (abs($Classified_ORFs->{$gene}->{$ORF}->{dist_aTIS}) <= 3 ) {$count_3_codon++}
				$count_Truncation++;
			} elsif ($Classified_ORFs->{$gene}->{$ORF}->{orf_type} eq "5' extension") {
				if (abs($Classified_ORFs->{$gene}->{$ORF}->{dist_aTIS}) <= 3) {$count_3_codon++}
				$count_Extension++;
			} elsif ($Classified_ORFs->{$gene}->{$ORF}->{orf_type} eq 'Intergenic') {
				$count_Intergenic++;
			} elsif ($Classified_ORFs->{$gene}->{$ORF}->{orf_type} eq "3' truncation") {
				$count_Internal++;
			} elsif ($Classified_ORFs->{$gene}->{$ORF}->{orf_type} eq 'Opposite strand') {
				$count_Reverse++;
			} elsif ($Classified_ORFs->{$gene}->{$ORF}->{orf_type} eq 'Out-of-frame') {
				$count_out_of_frame++;
			}else {
				$count_ncRNA++;
			}		
		}
	}

	print "Total number of ORFs Identified ", scalar(keys %$Classified_ORFs), "\n";
	print "number of close proximity variants $count_3_codon\n";
	print "Number of annotated ORFs $count_Annotated\n";
	print "Number of truncations ORFs $count_Truncation\n";
	print "Number of 5' extension ORFs $count_Extension\n";
	print "Number of Intergenic ORFs $count_Intergenic\n";
	print "Number of 3' truncation $count_Internal\n";
	print "Number of Reverse ORFs $count_Reverse\n";
	#print "Number of Out-of-frame ORFs $count_out_of_frame\n";
	print "Number of ncRNA ORFs $count_ncRNA\n";

} else {
	create_bed($ORFs_processed, $output_prefix);
	output_ORFs($ORFs_processed, $output_prefix);
	predicted_ORFs_pep($ORFs_processed,$output_prefix);
}



##################
##	SUBS
##################


sub ORF_classification {

	my $ORFs = $_[0];
	my $annotated = $_[1];

	foreach my $gene (sort keys %$ORFs) {

		my $ORF = (keys %{$ORFs->{$gene}})[0];
		my $start = $ORFs->{$gene}->{$ORF}->{start}; 
		my $stop = $ORFs->{$gene}->{$ORF}->{stop};
		my $strand = $ORFs->{$gene}->{$ORF}->{strand};
		my $region = $ORFs->{$gene}->{$ORF}->{region};

		if (exists($annotated->{$gene})) {
			if (exists $annotated->{$gene}->{$ORF}) {	# gene and ORF in annotation 
				if ($annotated->{$gene}->{$ORF}->{biotype} eq 'protein_coding') {
					$ORFs->{$gene}->{$ORF}->{orf_type} = 'Annotated';
					$ORFs->{$gene}->{$ORF}->{ref_anno} = $annotated->{$gene}->{$ORF}->{name};
					$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
				} elsif ($annotated->{$gene}->{$ORF}->{biotype} ne 'protein_coding')  {
					$ORFs->{$gene}->{$ORF}->{orf_type} = $annotated->{$gene}->{$ORF}->{biotype};
					$ORFs->{$gene}->{$ORF}->{ref_anno} = $annotated->{$gene}->{$ORF}->{name};
					$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
				}
			} else { # gene in anotation but not the ORF
				my $ORF_anno = (keys %{$annotated->{$gene}})[0];

				if (exists $annotated->{$gene}->{$ORF_anno} and $annotated->{$gene}->{$ORF_anno}->{biotype} eq 'protein_coding') {
					my $start_anno  = $annotated->{$gene}->{$ORF_anno}->{start}; 
					my $stop_anno   = $annotated->{$gene}->{$ORF_anno}->{stop};

					my $dist = ($strand eq '+') ? $start - $start_anno : $stop_anno - $stop;

					$ORFs->{$gene}->{$ORF}->{ref_anno} = $ORF_anno;
					$ORFs->{$gene}->{$ORF}->{dist_aTIS} = $dist;
					$ORFs->{$gene}->{$ORF}->{ref_anno} = $annotated->{$gene}->{$ORF_anno}->{name};

					if ($dist < 0) {
						$ORFs->{$gene}->{$ORF}->{orf_type} = "5' extension";
					} elsif ($dist > 0) {
						$ORFs->{$gene}->{$ORF}->{orf_type} = 'Truncation';
					} elsif ($dist == 0) {
						$ORFs->{$gene}->{$ORF}->{orf_type} = 'Annotated';
						$ORFs->{$gene}->{$ORF}->{ref_anno} = $annotated->{$gene}->{$ORF_anno}->{name};
						$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
					}
				} elsif (exists $annotated->{$gene}->{$ORF_anno} and $annotated->{$gene}->{$ORF_anno}->{biotype} ne 'protein_coding') {
					$ORFs->{$gene}->{$ORF}->{orf_type} = $annotated->{$gene}->{$ORF_anno}->{biotype};
					$ORFs->{$gene}->{$ORF}->{ref_anno} = $annotated->{$gene}->{$ORF_anno}->{name};
					$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
				} 
			} 
		} else {
			foreach my $gene_anno (sort keys %$annotated) {

				my $ORF_anno 	= (sort keys %{$annotated->{$gene_anno}})[0];
				my $strand_anno = $annotated->{$gene_anno}->{$ORF_anno}->{strand};
				my $region_anno = $annotated->{$gene_anno}->{$ORF_anno}->{region};
				my $start_anno  = $annotated->{$gene_anno}->{$ORF_anno}->{start}; 
				my $stop_anno   = $annotated->{$gene_anno}->{$ORF_anno}->{stop};
				my $biotype		= $annotated->{$gene_anno}->{$ORF_anno}->{biotype};
				my $name		= $annotated->{$gene_anno}->{$ORF_anno}->{name};

				next if ($region ne $region_anno);

				if ($start_anno == $start and $stop == $stop_anno and $biotype eq 'protein_coding')  {
					$ORFs->{$gene}->{$ORF}->{orf_type} = 'Annotated';
					$ORFs->{$gene}->{$ORF}->{ref_anno} = $name;
					$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
					last;
				} elsif ($start_anno <= $start and $stop <= $stop_anno and $biotype eq 'protein_coding') {
					if ($strand eq $strand_anno) {
						if ($start_anno == $start or $stop == $stop_anno ) {
							$ORFs->{$gene}->{$ORF}->{orf_type} = "3' truncation";
							$ORFs->{$gene}->{$ORF}->{ref_anno} = $name;
							$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
							last;
						} else {
							$ORFs->{$gene}->{$ORF}->{orf_type} = "Out-of-frame";
							$ORFs->{$gene}->{$ORF}->{ref_anno} = $name;
							$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
							last;
						}
					} else {
						$ORFs->{$gene}->{$ORF}->{orf_type} = "Opposite strand";
						$ORFs->{$gene}->{$ORF}->{ref_anno} = $name;
						$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
						last;
					}
				} elsif ($start_anno <= $start and $stop <= $stop_anno and $biotype ne 'protein_coding') {
					$ORFs->{$gene}->{$ORF}->{orf_type} = $biotype;
					$ORFs->{$gene}->{$ORF}->{ref_anno} = $name;
					$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
					last;
				}
			}
		}

		if (exists $ORFs->{$gene}->{$ORF} and defined $ORFs->{$gene}->{$ORF}->{orf_type}) {
			# do nothing
		} else {
			$ORFs->{$gene}->{$ORF}->{orf_type} = "Intergenic";
			$ORFs->{$gene}->{$ORF}->{ref_anno} = "NA";
			$ORFs->{$gene}->{$ORF}->{dist_aTIS} = 'NA';
			$ORFs->{$gene}->{$ORF}->{ref_anno} = 'NA';
		}
	}  # end of loop ove all ORFs

	return $ORFs;
}


sub ORFs_start_selection {

	my $ORFs = $_[0];

	my $ORFs_forward = {};	# hash to start selected forward strand ORfs
	my $ORFs_reverse = {};	# hash to start selected reverse strand ORfs

	# we group all ORFs by region => strand and stop
	my $forward_strand = {};
	my $reverse_strand = {};
	foreach my $gene (sort keys %$ORFs) {
		my $ORF = (keys %{$ORFs->{$gene}})[0];
		my $strand = $ORFs->{$gene}->{$ORF}->{strand};
		my $region = $ORFs->{$gene}->{$ORF}->{region};
		my $start = $ORFs->{$gene}->{$ORF}->{start};
		my $stop = $ORFs->{$gene}->{$ORF}->{stop};
		if ($strand eq '+') {
			$forward_strand->{$region}->{$stop} = $gene;
		} elsif ($strand eq '-') {
			$reverse_strand->{$region}->{$start} = $gene;
		}
	}

	# process forward strand
	foreach my $region (sort keys %$forward_strand) {
		my $strand = '+';
		my @all_stops = sort {$a<=>$b} keys %{$forward_strand->{$region}};	# sort gene families by their stop positions

		for (my $i=0; $i < scalar(@all_stops); $i++) {

			my $gene = $forward_strand->{$region}->{$all_stops[$i]};
			my $gene_orfs = get_all_ORFs($gene,$strand,$region);	
			my @gene_starts = sort {$b <=> $a} keys %$gene_orfs;	# for each gene family sort by start positions

			my $sel_idx = 0;
			for (my $j=1; $j <scalar(@gene_starts); $j=$j+1) {	
				# scan from the shortest and check if RPF reads upstream of shorter start
				# also take into account read accumulation within the SD region
				# If the start region overlaps an upstream gene then check peakedness using kurtosis measure

				my $start_sel = $gene_starts[$sel_idx];
				my $ORF_sel = $gene_orfs->{$start_sel};
				my $count_sel = $ORFs_RF->{$gene}->{$ORF_sel}->{ribo_count};

				my $start_cur = $gene_starts[$j];
				my $ORF_cur = $gene_orfs->{$start_cur};
				my $count_cur = $ORFs_RF->{$gene}->{$ORF_cur}->{ribo_count};
				my $avg_cur = $count_cur/$ORFs_RF->{$gene}->{$ORF_cur}->{len};

				if ($count_sel < $count_cur) {
					my ($flag_cur,$kurtosis_cur, $overlap_edge_cur) = check_overlap_and_kurtosis($start_cur,$j,$strand,$region,$forward_strand,\@all_stops);
					my ($flag_sel, $kurtosis_sel, $overlap_edge_sel) = check_overlap_and_kurtosis($start_sel,$sel_idx,$strand,$region,$forward_strand,\@all_stops);

					if ($flag_sel == 1 and $flag_cur == 1) {
						if ($kurtosis_cur > 0) {
							$sel_idx = $j;
						} else {
							if ($kurtosis_cur >= $kurtosis_sel) {$sel_idx = $j;}
						}
					} elsif ($flag_sel == 1 and $flag_cur == 0) {
						$sel_idx = $j;
					} elsif ($flag_sel == 0 and $flag_cur == 0) {
						my $SD_region_start = $ORFs_RF->{$gene}->{$ORF_sel}->{SD_pos} + $SDLEN;
						if ($start_cur - $SD_region_start > 0) {
							my ($count,$coverage,$rpkm) = RIBO_count($region,$strand,$SD_region_start-1,$start_cur);
							if ($count > $avg_cur) {$sel_idx = $j;}
						} else {
							$sel_idx = $j;
						}
					} elsif ($flag_sel == 0 and $flag_cur == 1) {
						if ($kurtosis_cur >= 0) {
							$sel_idx = $j;
						} else {
							my ($count,$coverage,$rpkm) = RIBO_count($region,$strand,$start_sel+3,$overlap_edge_sel);
							if ($count > $avg_cur) {$sel_idx = $j;}
						}
					}
				}
			}

			my $ORF_sel = $gene_orfs->{$gene_starts[$sel_idx]};
			$ORFs_reverse->{$gene}->{$ORF_sel} = $ORFs->{$gene}->{$ORF_sel};
		}
	}

	# process reverse strand
	foreach my $region (sort keys %$reverse_strand) {
		my $strand = '-';
		my @all_starts = sort {$a<=>$b} keys %{$reverse_strand->{$region}};

		for (my $i=scalar(@all_starts)-1; $i>=0; $i--) {

			my $gene = $reverse_strand->{$region}->{$all_starts[$i]};
			my $gene_orfs = get_all_ORFs($gene,$strand,$region);
			my @gene_stops = sort {$a <=> $b} keys %$gene_orfs;

			my $sel_idx = 0;
			for (my $j=1; $j <scalar(@gene_stops); $j++) {

				my $ORF_sel = $gene_orfs->{$gene_stops[$sel_idx]};
				my $start_sel = $gene_stops[$sel_idx];
				my $count_sel = $ORFs_RF->{$gene}->{$ORF_sel}->{ribo_count};

				my $ORF_cur = $gene_orfs->{$gene_stops[$j]};
				my $start_cur = $gene_stops[$j];
				my $count_cur = $ORFs_RF->{$gene}->{$ORF_cur}->{ribo_count};
				my $avg_cur = $count_cur/$ORFs_RF->{$gene}->{$ORF_cur}->{len};

				if ($count_sel < $count_cur) {
					my ($flag_cur,$kurtosis_cur, $overlap_edge_cur) = check_overlap_and_kurtosis($start_cur,$j,$strand,$region,$reverse_strand,\@all_starts);
					my ($flag_sel, $kurtosis_sel, $overlap_edge_sel) = check_overlap_and_kurtosis($start_sel,$sel_idx,$strand,$region,$reverse_strand,\@all_starts);
					if ($flag_sel == 1 and $flag_cur == 1) {
						if ($kurtosis_cur > 0) {
							$sel_idx = $j;
						} else {
							if ($kurtosis_cur >= $kurtosis_sel) {$sel_idx = $j;}
						}
					} elsif ($flag_sel == 1 and $flag_cur == 0) {
						$sel_idx = $j;
					} elsif ($flag_sel == 0 and $flag_cur == 0) {
						my $SD_region_start = $ORFs_RF->{$gene}->{$ORF_sel}->{SD_pos} + $SDLEN;
						if ($start_cur - $SD_region_start > 0) {
							my ($count,$coverage,$rpkm) = RIBO_count($region,$strand,$SD_region_start-1,$start_cur);
							if ($count > $avg_cur) {$sel_idx = $j;}
						} else {
							$sel_idx = $j;
						}
					} elsif ($flag_sel == 0 and $flag_cur == 1) {
						if ($kurtosis_cur >= 0) {
							$sel_idx = $j;
						} else {
							my ($count,$coverage,$rpkm) = RIBO_count($region,$strand,$start_sel+3,$overlap_edge_sel);
							if ($count > $avg_cur) {$sel_idx = $j;}
						}
					}
				}
			}

			my $ORF_sel = $gene_orfs->{$gene_stops[$sel_idx]};
			$ORFs_reverse->{$gene}->{$ORF_sel} = $ORFs->{$gene}->{$ORF_sel};
		}
	}

	my $ORFs_processed = {%$ORFs_forward, %$ORFs_reverse}; # combine both strand selections;

	return $ORFs_processed;

}


sub check_overlap_and_kurtosis {

	my $pos = $_[0];
	my $index = $_[1];
	my $strand = $_[2];
	my $region = $_[3];
	my $region_hash = $_[4];
	my @stops = @{$_[5]};
	
	# scan all downstreams ORFs to check for all overlapping ORFs
	my $flag = 0;
	my $overlap_edge;
	my $kurtosis = -100000; 

	if ($strand eq '+') {
		if ($index > 0) {
			my $escape = 0;
			for (my $i = $index - 1; $i >= 0; $i=$i-1) {
				if ($stops[$i] > $pos) {
					my $gene = $region_hash->{$region}->{$stops[$i]};
					foreach my $ORF (keys %{$ORFs_RF->{$gene}}) {
						my $start = $ORFs_RF->{$gene}->{$ORF}->{start};
						my $stop = $ORFs_RF->{$gene}->{$ORF}->{stop};
						if ($start < $pos and $pos < $stop) {
							$escape = 1;
							last;
						}
					}
					if ($escape == 1) {last}
				} else {
					last;
				}
			}
		}
	} else {
		if ($index < scalar(@stops) - 2) {
			my $escape = 0;
			for (my $i = $index + 1; $i < scalar(@stops); $i++) {
				if ($stops[$i] < $pos) {
					my $gene = $region_hash->{$region}->{$stops[$i]};
					foreach my $ORF (keys %{$ORFs_RF->{$gene}}) {
						my $start = $ORFs_RF->{$gene}->{$ORF}->{start};
						my $stop = $ORFs_RF->{$gene}->{$ORF}->{stop};
						if ($start < $pos and $pos < $stop) {
							$escape = 1;
							last;
						}
					}
					if ($escape == 1) {last}
				} else {
					last;
				}
			}
		}
	}

	# check kurtosis
	if ($flag == 1) {
		my @points = ();
		my $count_k = 0;
		my $WINDOW_UP = $WINDOW;
		my $start_k = ($strand eq '+') ? $pos - $WINDOW_UP: $pos - $WINDOW;
		my $stop_k = ($strand eq '+') ? $pos + $WINDOW: $pos + $WINDOW_UP;
		#for (my $t = $start_k; $t <= $stop_k + $WINDOW; $t=$t + 1) {
		for (my $t = $start_k; $t <= $stop_k + $WINDOW; $t=$t + 3) {
			if (defined $reads_table->{$region}->{$strand}->{$t}) {
				push @points, $reads_table->{$region}->{$strand}->{$t};
				$count_k += $reads_table->{$region}->{$strand}->{$t};
			} else {
				push @points, 0;
			}
		}
		if ($count_k > 0) {
			$kurtosis = kurtosis(@points);
		} 
	} else {
		$kurtosis = 1000;
	}

	return $flag, $kurtosis, $overlap_edge;
}


sub get_all_ORFs {

	my $gene = $_[0];
	my $strand = $_[1];
	my $region = $_[2];

	my $ORF_starts = {};
	if ($strand eq '+') {
		foreach my $ORF (keys %{$ORFs_RF->{$gene}}) {
			my $start = $ORFs_RF->{$gene}->{$ORF}->{start};
			$ORF_starts->{$start} = $ORF;
		}
	} else {
		foreach my $ORF (keys %{$ORFs_RF->{$gene}}) {
			my $stop = $ORFs_RF->{$gene}->{$ORF}->{stop};
			$ORF_starts->{$stop} = $ORF;
		}
	}

	return $ORF_starts;
}


sub kurtosis {

	my @nums = @_;

	my $sum = 0;
	foreach (@nums) { $sum += $_ }
	my $n = scalar(@nums);
	my $mean = $sum/$n;
	my $variance = 0;
	my $kurtosis = 0;
	foreach (@nums) {
		my $deviation = $_ - $mean;
		$variance += $deviation**2;
		$kurtosis += $deviation**4;
	}
	$variance /= ($n - 1);
	if ($variance) {
		$kurtosis = ($kurtosis/($n * $variance * $variance)) - 3.0;
	}

	return $kurtosis;
}


sub one_way_overlap {

	# calculate coverage and read density of the non-overlaping and overlapping regions
	# take a 30bps window (~ Ribosome length to account for any 5' of 3' mapping errors
	# drop ORF with non overlapping region less than MINCOV and/or MINREAD
	my $ORFs = $_[0];

	my $forward_strand = {};
	my $reverse_strand = {};
	foreach my $gene (sort keys %$ORFs) {
		my $ORF = (keys %{$ORFs->{$gene}})[0];
		my $strand = $ORFs->{$gene}->{$ORF}->{strand};
		my $region = $ORFs->{$gene}->{$ORF}->{region};
		my $start = $ORFs->{$gene}->{$ORF}->{start};
		my $stop = $ORFs->{$gene}->{$ORF}->{stop};
		if ($strand eq '+') {
			$forward_strand->{$region}->{$start}->{orf} = $ORF;
			$forward_strand->{$region}->{$start}->{gene} = $gene;
		} elsif ($strand eq '-') {
			$reverse_strand->{$region}->{$start}->{orf} = $ORF;
			$reverse_strand->{$region}->{$start}->{gene} = $gene;
		}
	}

	my $ORFs_to_delete = {};

	# foward strand
	foreach my $region (sort keys %$forward_strand) {
		my $strand = '+';
		my @all_stops = sort {$a<=>$b} keys %{$forward_strand->{$region}};	# sort gene families by their stop positions

		foreach my $start1 (sort {$a <=> $b} keys %{$forward_strand->{$region}}) {
			my $ORF1 = $forward_strand->{$region}->{$start1}->{orf};
			my $gene1 = $forward_strand->{$region}->{$start1}->{gene};
			my $stop1 = $ORFs->{$gene1}->{$ORF1}->{stop};

			next if (exists $ORFs_to_delete->{$ORF1});

			foreach my $start2 (sort {$a <=> $b} keys %{$forward_strand->{$region}}) {
				my $ORF2 = $forward_strand->{$region}->{$start2}->{orf};
				my $gene2 = $forward_strand->{$region}->{$start2}->{gene};
				my $stop2 = $ORFs->{$gene2}->{$ORF2}->{stop};

				next if (exists $ORFs_to_delete->{$ORF2});

				if ($start1 < $start2 and $start2 < $stop1 and $stop1 < $stop2) {
					my $length0 = $stop1 - $start2;	# length of overlapping region
					my $length1 = $start2 - $start1;
					my $length2 = $stop2 - $stop1;
					my $SD_start = $ORFs->{$gene2}->{$ORF2}->{SD_pos} - $SDLEN;
					my ($count0,$coverage0,$rpkm0) = RIBO_count($region,$strand,$start2,$stop1);
					my ($count1,$coverage1,$rpkm1) = RIBO_count($region,$strand,$start1,$SD_start-1);
					my ($count2,$coverage2,$rpkm2) = RIBO_count($region,$strand,$stop1+1,$stop2);

					if ($start1 < $SD_start) {
						if ($coverage1 < $MINCOV or $rpkm1 < $MINRPKM ) {
							if ($coverage2 >= $MINCOV and $rpkm2 >= $MINRPKM) {
								$ORFs_to_delete->{$ORF1} = 1; 
							} else {
								if ($length1 < $length0 and $length2 > $length0) {
									$ORFs_to_delete->{$ORF1} = 1;
								} elsif ($length2 < $length0 and $length1 > $length0) {
									$ORFs_to_delete->{$ORF2} = 1;
								}
							}
						} else {
							if ($coverage2 < $MINCOV or $rpkm2 < $MINRPKM) {
								$ORFs_to_delete->{$ORF2} = 1; 
							}
						}
					} else {
						# start region of ORF1 falls within SD region of ORF2
						if ($coverage2 < $MINCOV or $rpkm2 < $MINRPKM) {
							$ORFs_to_delete->{$ORF2} = 1; 
						} else {
							$ORFs_to_delete->{$ORF1} = 1; 
						}
					}
				}
			}
		}
	}

	# reverse strand
	foreach my $region (sort keys %$reverse_strand) {
		my $strand = '-';
		my @all_stops = sort {$a<=>$b} keys %{$reverse_strand->{$region}};	# sort gene families by their stop positions

		foreach my $start1 (sort {$a <=> $b} keys %{$reverse_strand->{$region}}) {
			my $ORF1 = $reverse_strand->{$region}->{$start1}->{orf};
			my $gene1 = $reverse_strand->{$region}->{$start1}->{gene};
			my $stop1 = $ORFs->{$gene1}->{$ORF1}->{stop};

			next if (exists $ORFs_to_delete->{$ORF1});

			foreach my $start2 (sort {$a <=> $b} keys %{$reverse_strand->{$region}}) {
				my $ORF2 = $reverse_strand->{$region}->{$start2}->{orf};
				my $gene2 = $reverse_strand->{$region}->{$start2}->{gene};
				my $stop2 = $ORFs->{$gene2}->{$ORF2}->{stop};

				next if (exists $ORFs_to_delete->{$ORF2});
				last if ($stop1 < $start2);

				if ($start1 < $start2 and $start2 < $stop1 and $stop1 < $stop2) {
					my $length0 = $stop1 - $start2;	# length of overlapping region
					my $length1 = $start2 - $start1;
					my $length2 = $stop2 - $stop1;
					my $SD_stop = $ORFs->{$gene1}->{$ORF1}->{SD_pos} + $SDLEN;
					my ($count0,$coverage0,$rpkm0) = RIBO_count($region,$strand,$start2-3,$stop1);
					my ($count1,$coverage1,$rpkm1) = RIBO_count($region,$strand,$start1,$start2-3);
					my ($count2,$coverage2,$rpkm2) = RIBO_count($region,$strand,$SD_stop+1,$stop2);
					if ($SD_stop < $stop2) {
						if ($coverage1 < $MINCOV or $rpkm1 < $MINRPKM) {
							if ($coverage2 >= $MINCOV and $rpkm2 >= $MINRPKM) {
								$ORFs_to_delete->{$ORF1} = 1; 
							} else {
								if ($length1 < $length0 and $length2 > $length0) {
									$ORFs_to_delete->{$ORF1} = 1;
								} elsif ($length2 < $length0 and $length1 > $length0) {
									$ORFs_to_delete->{$ORF2} = 1;
								}
							}
						} else {
							if ($coverage2 < $MINCOV or $rpkm2 < $MINRPKM) {
								$ORFs_to_delete->{$ORF2} = 1; 
							}
						}
					} else {
						# start region of ORF2 falls within SD region of ORF1
						if ($coverage1 < $MINCOV or $rpkm1 < $MINRPKM) {
							$ORFs_to_delete->{$ORF1} = 1; 
						} else {
							$ORFs_to_delete->{$ORF2} = 1; 
						}
					}
				}
			}
		}
	}

	my $selected_ORFs = {};
	foreach my $gene (keys %$ORFs) {
		foreach my $ORF (keys %{$ORFs->{$gene}}) {
			if (exists $ORFs_to_delete->{$ORF}) {
			} else {
				$selected_ORFs->{$gene}->{$ORF} = $ORFs->{$gene}->{$ORF};
			}
		}
	}

	return $selected_ORFs;
}


sub two_way_overlap {

	my $ORFs = $_[0];
	my $ORFs_to_delete = {};

	my $ORFs_by_pos = {};
	foreach my $gene (sort keys %$ORFs) {
		foreach my $ORF (sort keys %{$ORFs->{$gene}}) {
			my $region = $ORFs->{$gene}->{$ORF}->{region};
			my $start = $ORFs->{$gene}->{$ORF}->{start};
			my $strand = $ORFs->{$gene}->{$ORF}->{strand};
			$ORFs_by_pos->{$region}->{$strand}->{$start} = $ORF."#".$gene;
		}
	}

	my $two_count = 0;
	my $overlap = {};
	foreach my $region (keys %$ORFs_by_pos) {
		foreach my $strand (keys %{$ORFs_by_pos->{$region}}) {
			my @all_starts = sort {$a <=> $b} keys %{$ORFs_by_pos->{$region}->{$strand}};
			for (my $i = 1; $i < scalar(@all_starts) - 1; $i++) {
				last if ($i + 2 == scalar(@all_starts));

				# previous ORF
				my $start1 	= $all_starts[$i-1];
				my @orf_info1 = split('#', $ORFs_by_pos->{$region}->{$strand}->{$start1});
				my $ORF1  = $orf_info1[0];
				my $gene1 = $orf_info1[1];
				my $strand1 = $orf_info1[2];
				my ($stop1) = $ORF1 =~ /-(\d+)$/;
				$start1 = ($strand eq '+') ? $start1: $start1 - 3;
				$stop1 = ($strand eq '+') ? $stop1 + 3: $stop1;

				# Middle ORF
				my $start2 	= $all_starts[$i];
				my @orf_info2 = split('#', $ORFs_by_pos->{$region}->{$strand}->{$start2});
				my $ORF2  = $orf_info2[0];
				my $gene2 = $orf_info2[1];
				my ($stop2) = $ORF2 =~ /-(\d+)$/;
				$start2 = ($strand eq '+') ? $start2: $start2 - 3;
				$stop2 = ($strand eq '+') ? $stop2 + 3: $stop2;
			
				# preceeding ORF
				my $start3 	= $all_starts[$i+1];
				my @orf_info3 = split('#', $ORFs_by_pos->{$region}->{$strand}->{$start3});
				my $ORF3  = $orf_info3[0];
				my $gene3 = $orf_info3[1];
				my ($stop3) = $ORF3 =~ /-(\d+)$/;
				$start3 = ($strand eq '+') ? $start3: $start3 - 3;
				$stop3 = ($strand eq '+') ? $stop3 + 3: $stop3;

				# skip if orf before or after is in the deleted list
				next if (exists $ORFs_to_delete->{$ORF1} or exists $ORFs_to_delete->{$ORF3});
				if ($start1 < $start2 and $start2 < $stop1 and $stop1 < $stop2 and $start3 < $stop2 and $stop2 < $stop3) {
					$two_count++;
					if ($stop1 >= $start3 or abs($stop1 - $start3) < $REGION) {
						$ORFs_to_delete->{$ORF2} = 1;
					} else {
						# get mid region between two ORFs
						my $start_mid = $stop1 + $REGION;
						my $stop_mid = $start3 - $REGION;
						# check coverage and read density of entire mid point
						# if it is less than min cutoffs then drop ORF2else account for error in mapping by looking outside the ~30bps
						my ($count,$coverage,$rpkm)= (0,0,0);
						if ($stop_mid - $start_mid>0) {
							($count,$coverage,$rpkm) = RIBO_count($region,$strand,$start_mid,$stop_mid);
							if  ($rpkm < $MINRPKM or $coverage < $MINCOV) {
								$ORFs_to_delete->{$ORF2} = 1;
							}
						} else {
							$ORFs_to_delete->{$ORF2} = 1;
						}
					}
				}
			}
		}
	}

	my $selected_ORFs = {};
	foreach my $gene (keys %$ORFs) {
		foreach my $ORF (keys %{$ORFs->{$gene}}) {
			unless (exists $ORFs_to_delete->{$ORF}) {
				$selected_ORFs->{$gene}->{$ORF} = $ORFs->{$gene}->{$ORF};
			}
		}
	}

	return $selected_ORFs;
}


sub RIBO_count {

	my $chr 	  = $_[0];
	my $strand 	  = $_[1];
	my $start 	  = $_[2];
	my $stop 	  = $_[3];

	my $count = 0;
	my $coverage = 0;
	my $rpkm = 0;
	for ( my $pos = $start; $pos <= $stop; $pos++) {
		if ($reads_table->{$chr}->{$strand}->{$pos}) {
			$count += $reads_table->{$chr}->{$strand}->{$pos};
			$coverage++;
		}
	}
	my $length = $stop - $start + 1;
	if ($length <= 0) {$count = 0; $coverage=0; $length=1} else {$coverage = $coverage/$length;}
	$rpkm = ($count*1000000000)/($length*$mapped_total);

	return $count,$coverage,$rpkm;

}


sub internal_out_of_frame {
	# find all internal out of frame ORFs
	my $ORFs = $_[0];

	my $ORFs_to_delete = {};
	my $count = 0;
	foreach my $gene1 (sort keys %$ORFs) {
		my $ORF1 = (keys %{$ORFs->{$gene1}})[0];
		next if (exists $ORFs_to_delete->{$ORF1});

		my $strand1 = $ORFs->{$gene1}->{$ORF1}->{strand};
		my $start1 = $ORFs->{$gene1}->{$ORF1}->{start};
		my $stop1 = $ORFs->{$gene1}->{$ORF1}->{stop};
		my $region1 = $ORFs->{$gene1}->{$ORF1}->{region};

		foreach my $gene2 (sort keys %$ORFs) {
			next if ($gene1 eq $gene2);
			my $ORF2 = (keys %{$ORFs->{$gene2}})[0];
			next if ($ORF1 eq $ORF2);
			next if (exists $ORFs_to_delete->{$ORF2});

			my $strand2 = $ORFs->{$gene2}->{$ORF2}->{strand};
			my $start2 = $ORFs->{$gene2}->{$ORF2}->{start};
			my $stop2 = $ORFs->{$gene2}->{$ORF2}->{stop};
			my $region2 = $ORFs->{$gene2}->{$ORF2}->{region};

			if ($strand1 eq $strand2 and $region1 eq $region2) {	# consider only ORFs on same strand and region
				if ($start1 < $start2 and $stop2 < $stop1) {
					$ORFs_to_delete->{$ORF2} = 1;
					$count++;
				}
			}
		}
	}

	my $selected_ORFs = {};
	foreach my $gene (keys %$ORFs) {
		foreach my $ORF (keys %{$ORFs->{$gene}}) {
			next if (exists $ORFs_to_delete->{$ORF});
			$selected_ORFs->{$gene}->{$ORF} = $ORFs->{$gene}->{$ORF};
		}
	}

	return $selected_ORFs;
}


sub Accumulation_proportion {

	my $start = $_[0];
	my $stop = $_[1];
	my $strand = $_[2];
	my $region = $_[3];

	my $length = $stop - $start + 1;
	my $tmp_OFFSET_START = $OFFSET_START;

	if ($length < $OFFSET_START) {$tmp_OFFSET_START = ceil(0.7*$length)}
	my $start_start = ($strand eq '+') ? $start - 3: $stop - $tmp_OFFSET_START;
	my $start_stop = ($strand eq '+') ? $start + $tmp_OFFSET_START: $stop+3;
	my $rest_start = ($strand eq '+') ? $start + $tmp_OFFSET_START + 1: $start;
	my $rest_stop = ($strand eq '+') ? $stop : $stop - $tmp_OFFSET_START - 1;

	# average reads
	my $average = 0;
	for ( my $pos = $start; $pos <= $stop; $pos++) {
		if ($reads_table->{$region}->{$strand}->{$pos}) {
			$average += $reads_table->{$region}->{$strand}->{$pos};
		} 
	}
	$average = $average/($stop - $start + 1);

	my $count_rest = 0;
	my $count_start = 0;
	for ( my $pos = $start; $pos <= $stop; $pos++) {
		if ($start_start <= $pos and $pos <= $start_stop) { 
			if ($reads_table->{$region}->{$strand}->{$pos}) {
				$count_start += $reads_table->{$region}->{$strand}->{$pos};
			} else {
				$count_start += $average;
			}
		}

		if ($rest_start <= $pos and $pos <= $rest_stop) { 
			if ($reads_table->{$region}->{$strand}->{$pos}) {
				$count_rest += $reads_table->{$region}->{$strand}->{$pos};
			} else {
				$count_rest += $average;
			}
		}
	}

	my $average_start = $count_start/(abs($start_stop - $start_start) + 1);
	my $average_rest = $count_rest/(abs($rest_stop - $rest_start) + 1);

	if ($average_rest == 0) {
		return 0, $average_start, $average_rest;
	} else {
		return $average_start/$average_rest, $average_start, $average_rest;
	}
}


sub output_ORFs_classified {

	my $ORFs = $_[0];
	my $prefix = $_[1];

	open (F, ">".$predicted_ORFs) or die  "Error creating file: $predicted_ORFs\n";
	print F "ORF_locus\tstrand\tlength\tstart_codon\tribo_count\tribo_rpkm\tribo_coverage\tSD_score\tSD_pos\tprob\tORF_type\tReference\tDistance_from_aTIS\n";
	foreach my $gene (sort keys %$ORFs) {
		foreach my $ORF (sort keys %{$ORFs->{$gene}}) {
			my $strand = $ORFs->{$gene}->{$ORF}->{strand};
			my $length = $ORFs->{$gene}->{$ORF}->{len};
			my $start_codon = $ORFs->{$gene}->{$ORF}->{start_codon}; 
			my $ribo_count = $ORFs->{$gene}->{$ORF}->{ribo_count}; 
			my $ribo_rpkm = $ORFs->{$gene}->{$ORF}->{ribo_rpkm};
			my $ribo_coverage = $ORFs->{$gene}->{$ORF}->{coverage};
			my $SD_score = $ORFs->{$gene}->{$ORF}->{SD_score}; 
			my $SD_pos = $ORFs->{$gene}->{$ORF}->{SD_pos}; 
			my $prob = $ORFs->{$gene}->{$ORF}->{prob}; 
			my $region = $ORFs->{$gene}->{$ORF}->{region};
			my $orf_type = $ORFs->{$gene}->{$ORF}->{orf_type};
			my $ref_anno = $ORFs->{$gene}->{$ORF}->{ref_anno};
			my $dist_aTIS = $ORFs->{$gene}->{$ORF}->{dist_aTIS};

			print F "$ORF\t$strand\t$length\t$start_codon\t$ribo_count\t$ribo_rpkm\t$ribo_coverage\t$SD_score\t$SD_pos\t$prob\t$orf_type\t$ref_anno\t$dist_aTIS\n";
		}
	}
	close F;
}

sub output_ORFs {

	my $ORFs = $_[0];
	my $prefix = $_[1];

	open (F, ">".$predicted_ORFs) or die  "Error creating file: $predicted_ORFs\n";
	print F "ORF_locus\tstrand\tlength\tstart_codon\tribo_count\tribo_rpkm\tribo_coverage\tSD_score\tSD_pos\tprob\n";
	foreach my $gene (sort keys %$ORFs) {
		foreach my $ORF (sort keys %{$ORFs->{$gene}}) {
			my $strand = $ORFs->{$gene}->{$ORF}->{strand};
			my $length = $ORFs->{$gene}->{$ORF}->{len};
			my $start_codon = $ORFs->{$gene}->{$ORF}->{start_codon}; 
			my $ribo_count = $ORFs->{$gene}->{$ORF}->{ribo_count}; 
			my $ribo_rpkm = $ORFs->{$gene}->{$ORF}->{ribo_rpkm};
			my $ribo_coverage = $ORFs->{$gene}->{$ORF}->{coverage};
			my $SD_score = $ORFs->{$gene}->{$ORF}->{SD_score}; 
			my $SD_pos = $ORFs->{$gene}->{$ORF}->{SD_pos}; 
			my $prob = $ORFs->{$gene}->{$ORF}->{prob}; 
			my $region = $ORFs->{$gene}->{$ORF}->{region};

			print F "$ORF\t$strand\t$length\t$start_codon\t$ribo_count\t$ribo_rpkm\t$ribo_coverage\t$SD_score\t$SD_pos\t$prob\n";
		}
	}
	close F;
}


sub create_bed {

	my $ORFs = $_[0];
	my $name = $_[1];

	my ($a,$b) = $name =~ /(.*)\/(.*)/;

	open (F, ">".$predicted_ORFs_bed) or die  "Error creating file: $predicted_ORFs_bed\n";
	print F "track type=bed name=\"$b\" description=\"\" visibility=1 color=\"7,7,255\" priority=20\n";
	foreach my $gene (sort keys %$ORFs) {
		my $ORF = (keys %{$ORFs->{$gene}})[0];
		my $region = $ORFs->{$gene}->{$ORF}->{region};
		my $start = $ORFs->{$gene}->{$ORF}->{start}; 
		my $stop = $ORFs->{$gene}->{$ORF}->{stop};
		my $strand = $ORFs->{$gene}->{$ORF}->{strand};
		my $orftype = $ORFs->{$gene}->{$ORF}->{orf_type};
		my $bed = "";
		if ($strand eq '+') {
			$start -= 1;
			my $thickStop = $stop;
			$stop = $stop + 3;
			$bed = "$region\t$start\t$stop\t$ORF\t1\t$strand\t$start\t$thickStop";
		} else {
			$start-=1;
			my $thickStart = $start;
			$start = $start - 3;
			$bed = "$region\t$start\t$stop\t$ORF\t1\t$strand\t$thickStart\t$stop";
		}
		print F "$bed\n";
	}  
	close F;
}


sub create_bed_classified {

	my $ORFs = $_[0];
	my $name = $_[1];

	my ($a,$b) = $name =~ /(.*)\/(.*)/;

	open (F, ">".$predicted_ORFs_bed) or die  "Error creating file: $predicted_ORFs_bed\n";
	print F "track type=bed name=\"$b\" description=\"\" visibility=1 color=\"\" priority=20\n";
	foreach my $gene (sort keys %$ORFs) {

		my $ORF = (keys %{$ORFs->{$gene}})[0];
		my $region = $ORFs->{$gene}->{$ORF}->{region};
		my $start = $ORFs->{$gene}->{$ORF}->{start}; 
		my $stop = $ORFs->{$gene}->{$ORF}->{stop};
		my $strand = $ORFs->{$gene}->{$ORF}->{strand};
		my $orftype = $ORFs->{$gene}->{$ORF}->{orf_type};

		my $itemRgb = "";
		if ($orftype eq "Annotated") {
			$itemRgb = "7,7,255";
		} elsif ($orftype eq "Extension") {
			$itemRgb = "138,6,255";
		} elsif ($orftype eq "Truncation") {
			$itemRgb = "4,214,255";
		} else {
			$itemRgb = "149,9,7";
		}

		my $bed = "";
		if ($strand eq '+') {
			$start -= 1;
			my $thickStop = $stop;
			$stop = $stop + 3;
			$bed = "$region\t$start\t$stop\t$ORF\t1\t$strand\t$start\t$thickStop\t$itemRgb";
		} else {
			$start-=1;
			my $thickStart = $start;
			$start = $start - 3;
			$bed = "$region\t$start\t$stop\t$ORF\t1\t$strand\t$thickStart\t$stop\t$itemRgb";
		}

		print F "$bed\n";
	}  # end of loop ove all ORFs
	close F;

	return $ORFs;

}


sub predicted_ORFs_pep {
	
	my $ORFs = $_[0];
	my $name = $_[1];

	open (F, ">".$predicted_ORFs_fasta) or die  "Error creating file: $predicted_ORFs_fasta\n";
	foreach my $gene (sort keys %$ORFs) {
		foreach my $ORF (sort keys %{$ORFs->{$gene}}) {
			my $strand = $ORFs->{$gene}->{$ORF}->{strand};
			my $region = $ORFs->{$gene}->{$ORF}->{region};
			my $start  = $ORFs->{$gene}->{$ORF}->{start}; 
			my $stop   = $ORFs->{$gene}->{$ORF}->{stop};
			my $length = $ORFs->{$gene}->{$ORF}->{len};
			my $start_codon = $ORFs->{$gene}->{$ORF}->{start_codon};

			my ($aa_seq,$dna_seq) = translate($start,$stop,$strand,$region,$ORF);
			$ORFs->{$gene}->{$ORF}->{aa_seq} = $aa_seq;

			$aa_seq =~ s/^./M/;
			print F ">generic|$ORF|start codon:$start_codon strand:$strand length:$length\n$aa_seq\n";
		}
	}
	close F;
}


sub translate {

	my $start = $_[0];
	my $stop = $_[1];
	my $strand = $_[2];
	my $region = $_[3];
	my $ORF = $_[4];

	my $length = $stop - $start + 1;
	my $dna_seq = ($strand eq '+') ? substr($genomes->{$region}, $start - 1,$length) : revdnacomp(substr($genomes->{$region}, $start-1,$length));

	my $aa_seq = "";
	for (my $i = 0; $i <= (length($dna_seq) - 3); $i = $i + 3) {
		my $codon = uc(substr($dna_seq, $i, 3));
		last if ($codon eq '*');
		$aa_seq = $aa_seq.$translationHash{$codon};
	}

	return $aa_seq, $dna_seq;

}


sub read_file {

	my $file = $_[0];
	my $ORFs = {};

	open (F, $file) or die  "Error reading file: $file";
	while (<F>) {
		chomp $_;
		if (/^orf_id/) {next;}
		my @line 			= (split '\t', $_);
		my $ORF 			= $line[0];
		my $gene			= $line[1];
		my $strand 			= $line[2];
		my $length 			= $line[3];
		my $start_codon 	= $line[4];
		my $ribo_count		= $line[5];
		my $ribo_rpkm		= $line[6];
		my $ribo_coverage	= $line[7];
		my $SD_score 		= $line[8];
		my $SD_pos			= $line[9];
		my $prob 			= $line[10];

		($ORFs->{$gene}->{$ORF}->{start}) 	= $ORF =~ /:(\d+)-/;
		($ORFs->{$gene}->{$ORF}->{stop}) 	= $ORF =~ /-(\d+)$/;
		my ($region) = $ORF =~ /^(.*):/;

		$ORFs->{$gene}->{$ORF}->{region} 	= $region;
		$ORFs->{$gene}->{$ORF}->{strand} 	= $strand;
		$ORFs->{$gene}->{$ORF}->{len} 		= $length;
		$ORFs->{$gene}->{$ORF}->{start_codon} = $start_codon;
		$ORFs->{$gene}->{$ORF}->{ribo_rpkm} = $ribo_rpkm;
		$ORFs->{$gene}->{$ORF}->{ribo_count}= $ribo_count;
		$ORFs->{$gene}->{$ORF}->{coverage} 	= $ribo_coverage;
		$ORFs->{$gene}->{$ORF}->{SD_score} 	= $SD_score;
		$ORFs->{$gene}->{$ORF}->{SD_pos} 	= $SD_pos;
		$ORFs->{$gene}->{$ORF}->{prob} 		= $prob;
	}
	close F; 

	return $ORFs;
}


sub annotation {

	my $file = $_[0];
	my $transcripts = {};
	my $annotated = {};

	open (F, $file) or die  "Error reading file: $file";
	while (<F>) {
		next if (/^#/);
		chomp $_;

		my @line = (split '\t', $_);
		my $region = $line[0];
		my $feature = $line[2];
		my $start = $line[3];
		my $stop = $line[4];
		my $strand = $line[6];

		if ($feature eq 'transcript') {

			my ($biotype) = $_ =~ /gene_biotype."?([^";]+)"?/;
			my ($gene_name) = $_ =~ /gene_name."?([^";]+)"?/;
			my ($gene_id) = $_ =~ /gene_id."?([^";]+)"?/;
			my ($tr_id) = $_ =~ /transcript_id."?([^";]+)"?/;

			my $name = ($gene_name) ? $gene_name : $gene_id;
			if ($biotype eq 'protein_coding') {
				if ($strand eq '+') {$stop = $stop - 3;} 
				else {$start = $start + 3;}
			}
		 	
			if (exists $annotated->{$tr_id}) {
				if ($start < $annotated->{$tr_id}->{start}) {
		 			$annotated->{$tr_id}->{start} = $start;
				}

				if ($stop > $annotated->{$tr_id}->{start}) {
		 			$annotated->{$tr_id}->{stop} = $stop;
				}
			} else {
		 		$annotated->{$tr_id}->{start} = $start;
		 		$annotated->{$tr_id}->{stop} = $stop;
		 		$annotated->{$tr_id}->{biotype} = $biotype;
		 		$annotated->{$tr_id}->{region} = $region;
		 		$annotated->{$tr_id}->{strand} = $strand;
		 		$annotated->{$tr_id}->{name} = $name;
			}
		}
	}
	close F;

	my $seen = {};

	foreach my $tr_id (keys %$annotated) {

		my $strand = $annotated->{$tr_id}->{strand};
		my $region = $annotated->{$tr_id}->{region};
		my $start = $annotated->{$tr_id}->{start};
		my $stop = $annotated->{$tr_id}->{stop};

		my $gene = ($strand eq "+") ? $region.":FAM+".$stop: $region.":FAM-".$start;
		my $chr_id = $region.":".$start."-".$stop;

	 	$transcripts->{$gene}->{$chr_id}->{biotype} = $annotated->{$tr_id}->{biotype};
	 	$transcripts->{$gene}->{$chr_id}->{region} = $region;
	 	$transcripts->{$gene}->{$chr_id}->{start} = $annotated->{$tr_id}->{start};
	 	$transcripts->{$gene}->{$chr_id}->{stop} = $annotated->{$tr_id}->{stop};
	 	$transcripts->{$gene}->{$chr_id}->{strand} = $strand;
	 	$transcripts->{$gene}->{$chr_id}->{name} = $annotated->{$tr_id}->{name};
	}

	return $transcripts;
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
	
	return $read,$mapped_total;
}


sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}


