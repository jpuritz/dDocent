#!/usr/bin/perl
use warnings;
use strict;

# Thanks to Simon Andrews for this code
# Merge together two FastQ files
# Usage is merge_fastq.pl [read1 file] [read2 file] [outfile]


my ($in1,$in2,$out) = @ARGV;

die "Usage is merge_fastq.pl [read1 file] [read2 file] [outfile]\n" unless ($out);

open (IN1,$in1) or die "Can't open $in1: $!";
open (IN2,$in2) or die "Can't open $in2: $!";
open (OUT,'>',$out) or die "Can't write to $out: $!";

my $count;
while (1) {
  ++$count;
  my $line1 = <IN1>;
  my $line2 = <IN2>;

  last unless (defined $line1 and defined $line2);

  if ($count % 2) {
    print OUT $line1;
  }
  else {
    chomp $line1;
    print OUT $line1,("N"x10),$line2;
  }

}
