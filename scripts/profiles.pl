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

my ($positive_set,$occupancyFile,$MINREAD,$work_dir,$script_dir) = @ARGV;

my $MINCOV = 0.1;
my $UPSTREAM = 10;
my $DOWNSTREAM = 300;

# read RPF into meory
my ($RPF,$mapped_total) = get_read_table($occupancyFile);

# read positive set
my $ORFs = read_file($positive_set);

my $start_file = $work_dir."tmp/start_profile.txt";
my $stop_file = $work_dir."tmp/stop_profile.txt";

my $tr_count = profile_data($ORFs,$RPF,$start_file,$stop_file);
system("Rscript ".$script_dir."/plot_profile.R $start_file $stop_file $work_dir $tr_count");



##################
##	SUBS
##################

sub profile_data {

	my $ORFs = $_[0];
	my $RPF = $_[1];
	my $start_file = $_[2];
	my $stop_file = $_[3];

	my @start_profile = (0)x($DOWNSTREAM + $UPSTREAM+1);
	my @stop_profile = (0)x($DOWNSTREAM + $UPSTREAM+1);
	my @start_count = (1)x($DOWNSTREAM + $UPSTREAM+1);
	my @stop_count = (1)x($DOWNSTREAM + $UPSTREAM+1);

	my $tr_count = 0;
	foreach my $tr (keys %$ORFs) {
		my $region = $ORFs->{$tr}->{region};
		my $strand = $ORFs->{$tr}->{strand};
		my $start = $ORFs->{$tr}->{start};
		my $stop = $ORFs->{$tr}->{stop};

		my $read_count = 0;
        my $coverage = 0;
		for (my $p = $start; $p <= $stop; $p++) {
			if ($RPF->{$region}->{$strand}->{$p}) {
				$read_count += $RPF->{$region}->{$strand}->{$p};
                $coverage++;
			} 
		}

		my $length = $stop-$start+1;
		#next if ($read_count < $MINREAD);
        next if (($coverage/$length) < $MINCOV);

		my $DOWNSTREAM_tmp = $DOWNSTREAM;
		if ($length < $DOWNSTREAM) {$DOWNSTREAM_tmp = $length}


		$tr_count++;

		# Start profile
		if ($strand eq '+') {
			my $t = 0;
			for (my $i = $start - $UPSTREAM; $i <= $start + $DOWNSTREAM_tmp; $i++) {
				if ($RPF->{$region}->{$strand}->{$i}) {
					$start_profile[$t] = $start_profile[$t] + $RPF->{$region}->{$strand}->{$i}/$read_count;
					$start_count[$t] = $start_count[$t] + 1;
				}
				$t++;
			}
		} else {
			my $t = 0;
			for (my $i = $stop + $UPSTREAM; $i >= $stop - $DOWNSTREAM_tmp; $i--) {
				if ($RPF->{$region}->{$strand}->{$i}) {
					$start_profile[$t] = $start_profile[$t] + $RPF->{$region}->{$strand}->{$i}/$read_count;
					$start_count[$t] = $start_count[$t] + 1;
				}
				$t++;
			}
		}

		# Stop profile
		if ($strand eq '+') {
			my $t = 0;
			for (my $i = $stop - $DOWNSTREAM_tmp; $i <= $stop + $UPSTREAM; $i++) {
				if ($RPF->{$region}->{$strand}->{$i}) {
					$stop_profile[$t] = $stop_profile[$t] + $RPF->{$region}->{$strand}->{$i}/$read_count;
					$stop_count[$t] = $stop_count[$t] + 1;
				}
				$t++;
			}
		} else {
			my $t = 0;
			for (my $i = $start + $DOWNSTREAM_tmp; $i >= $start - $UPSTREAM; $i--) {
				if ($RPF->{$region}->{$strand}->{$i}) {
					$stop_profile[$t] = $stop_profile[$t] + $RPF->{$region}->{$strand}->{$i}/$read_count;
					$stop_count[$t] = $stop_count[$t] + 1;
				}
				$t++;
			}
		}
	}

	open F1, ">".$start_file or die $!;
	print F1 "position\tAverage_read\n";
	my $p = -$UPSTREAM;
	for (my $i = 0; $i < scalar(@start_profile); $i++) {
		my $avg_reads = $start_profile[$i]/$start_count[$i];
		print F1 "$p\t$avg_reads\n";
		$p++;
	}
	close F1;

	open F2, ">".$stop_file or die $!;
	print F2 "position\tAverage_read\n";
	$p = -$DOWNSTREAM;
	foreach (my $i = 0; $i < scalar(@stop_profile); $i++) {
		my $avg_reads = $stop_profile[$i]/$stop_count[$i];
		print F2 "$p\t$avg_reads\n";
		$p++;
	}
	close F2;

	return $tr_count;
}


sub read_file {

	my $file = $_[0];
	my $ORFs = {};

	open (F, $file) or die  "Error reading file: $file";
	while (<F>) {
		chomp $_;
		next if (/^orf_id/);
		my @line = (split '\t', $_);
		my $ORF = $line[0];
		my $strand = $line[1];

		($ORFs->{$ORF}->{start}) = $ORF =~ /:(\d+)-/;
		($ORFs->{$ORF}->{stop}) = $ORF =~ /-(\d+)$/;
		($ORFs->{$ORF}->{region}) = $ORF =~ /^(.*):/;
		$ORFs->{$ORF}->{strand} = $strand;
	}
	close F; 

	return $ORFs;
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
			my @line = split ',|\t', $_;
			my $region 	= $line[0];
			my $start 	= $line[1];
			my $strand 	= $line[2];
			$read->{$region}->{$strand}->{$start} = $line[3];
		}
	}
	
	return $read,$mapped_total;
}

