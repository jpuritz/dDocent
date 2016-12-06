#!/usr/bin/perl -w

#Script written by Chris Hollenbeck.  Contact him at chollenbeck07@neo.tamu.edu

use strict;
use Getopt::Long;
use Pod::Usage;

pod2usage(-verbose => 1) if @ARGV == 0;

my $vcffile = '';
my $outfile = 'new.hwe';
my $popmap = '';
my $hwe = 0.001;
my $cutoff = 0.25;

GetOptions(	'vcf|v=s' => \$vcffile,
			'out|o=s' => \$outfile,
			'popmap|p=s' => \$popmap,
			'hwe|h=s' => \$hwe,
			'cutoff|c=s' => \$cutoff,
			);

unless ($vcffile) {
	print "\nNeed to specify a VCF file (-v) for input\n\n";
	pod2usage(-verbose => 1);
}

unless ($popmap) {
	print "\nNeed to specify a population map (-p) for input\n\n";
	pod2usage(-verbose => 1);
}


open(POP, "<", $popmap) or die $!;
my %pops;
while(<POP>) {
	next if $_ =~ /^\s/;
	chomp;
	my ($sample, $pop) = split;
	$pops{$pop} = [] unless $pops{$pop};
	push @{$pops{$pop}}, $sample;
}
close POP;

my %exclude_count;
foreach my $pop (sort keys %pops) {
	open(INDO, ">", $pop . '.inds') or die $!;
	foreach my $ind (@{$pops{$pop}}) {
		print INDO $ind, "\n";
	}
	close INDO;

	my $indfile = $pop . '.inds';

	print "Processing population: $pop (" , scalar(@{$pops{$pop}}) , " inds)", "\n";

	my $ouput = `vcftools --vcf $vcffile --keep $indfile --hardy --out $pop 2>&1`;

	open(HWEI, "<", $pop . '.hwe') or die $!;

	<HWEI>; 	# process the header
	while(<HWEI>) {
		last if $_ =~ /^\s/;
		chomp;
		my ($locus, $pos, $obs, $exp, $chisq, $pvalue, @rest) = split;
		$exclude_count{"$locus-$pos"}++ if $pvalue < $hwe;

	}
	close HWEI;

}

open(HWEO, ">", 'exclude.hwe') or die $!;
foreach my $snp (keys %exclude_count) {
	my ($locus, $site) = split('-', $snp);
	if ($exclude_count{$snp} / scalar(keys %pops) > $cutoff) {
		print HWEO join("\t", $locus, $site), "\n";
	}
}
close HWEO;

my $output = `vcftools --vcf $vcffile --exclude-positions exclude.hwe --recode --recode-INFO-all --out $outfile 2>&1`;
my $filt_output = `vcftools --vcf $vcffile --positions exclude.hwe --hardy --out filtered 2>&1`;

print "Outputting results of HWE test for filtered loci to 'filtered.hwe'\n";

my $kept;
my $total;
if ($output =~ /kept (\d+) out of a possible (\d+) Sites/) {
	$kept = $1;
	$total = $2;
}

print "Kept $kept of a possible $total loci (filtered " , $total - $kept , ' loci)', "\n";

__END__

=head1 NAME

filter_hwe_by_pop.pl

=head1 SYNOPSIS

filter_hwe_by_pop.pl -v <vcffile> -p <popmap> [options]

Options:
     -v     <vcffile>	input vcf file
     -p		<popmap>	tab-separated file of samples and population designations
	 -h		[hwe]	minimum Hardy-Weinberg p-value cutoff for SNPs
	 -c		[cutoff]	proportion of all populations that a locus can be below HWE cutoff without being filtered
     -o		[out]	name of outfile


=head1 OPTIONS

=over 8

=item B<-v, --vcffile>

VCF input file

=item B<-p, --popmap>

File with names of individuals and population designations, one per line

=item B<-h, --hwe>

Minimum cutoff for Hardy-Weinberg p-value (for test as implemented in vcftools) [Default: 0.001]

=item B<-c, --cutoff>

Proportion of all populations that a locus can be below HWE cutoff without being filtered. For example, choosing 0.5 will filter SNPs that are below the p-value threshold in 50% or more of the populations. [Default: 0.25]

=item B<-o, --out>

Name of outfile, by vcftools conventions (will be named X.recode.vcf)

=back

=head1 DESCRIPTION

B<filter_hwe_by_pop.pl> is a Perl wrapper for vcftools, designed to run tests for HWE on multiple populations and exclude loci that fall beneath a given threshold from the entire dataset

=cut
