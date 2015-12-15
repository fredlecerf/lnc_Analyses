#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# Name: lnc_01_parseBlastResults.pl
#-------------------------------------------------------------------------------
# Desc: Analyze BLAST results (XML output) from a simple command line
#
# The original version is FROM ACG_04_parseBlastResults
#
# VERSION HISTORY:
# 1.0 : initial release
# 1.01: bug fixes where coverage could > 100%. In initial version, the coverage
# was computed with hsp_length = $hsp->length('total'). This method includes the gaps
# to return the length of the HSP. Instead, the keyword 'hit' is now in place for
# version 1.01
# 1.02: adding query coverage and expect as input parameters
#-------------------------------------------------------------------------------
# Author: Frédéric lecerf on sept. 2013
#-------------------------------------------------------------------------------
use strict;
use Bio::SearchIO; 
use Getopt::Long;


#-------------------------------------------------------------------------------
# CONSTANT value definition
#-------------------------------------------------------------------------------
$|=1;
my $soft = "acg_04_parseBlastResults.pl";
my $VERSION = "1.02";
my $year = "2013";

# Benchmark 
my $begin_time = times();
my $verbose = undef;
my $inputFile = undef;
my $input_cov = undef;
my $input_exp = undef;

my @Options = @ARGV;

my %Options = (
        'input=s'	 	=> \$inputFile,
		'verbose'		=> \$verbose,
		'coverage=s'	=> \$input_cov,
		'expect=s'		=> \$input_exp,
	      );

my %OptionsHelp = (
		    'input    ' => '-i [file], XML blast output',
			'verbose  ' => '-v, verbose mode (optional)',
			'coverage ' => '-c [nb], query coverage threshold value (in %, optional, default: 70%)',
			'expect   ' => '-e [nb], add expect threshold (optional, default 10-5)',
		    );

my $param = join(' ',@Options);
#-------------------------------------------------------------------------------
sub usage ( $ ) {
	my ($msg)=@_;
	print STDERR "Error: $msg\n";
	print "--------------------------------------------------------------------------------\n";
	print "$soft $VERSION ($year)\n";
	print "--------------------------------------------------------------------------------\n";
	print "Desc: Using XML Blast output, parse and filter in a comprehensive format\n";
	print "for further analyses\n";
	print "--------------------------------------------------------------------------------\n";
	print STDERR "Usage: $0 [options]\n";
	print STDERR "see README for specific configuration parameters\n";
	print STDERR "Options:\n";
	map {printf STDERR "\t$_: %s\n",$OptionsHelp{$_};} keys %OptionsHelp; 
	print STDERR "Please cite: Lecerf F., $year\n";
	exit(1);
}


#-------------------------------------------------------------------------------
# Check parameters & read data file
#-------------------------------------------------------------------------------
GetOptions(%Options);
defined $inputFile			|| &usage('XML blast output required (-i)');

print "--------------------------------------------------------------------------------\n";
print "$soft $VERSION ($year)\n";
print "--------------------------------------------------------------------------------\n";

print "Loading file: $inputFile\n";
my $in = new Bio::SearchIO(-format => 'blastxml', 
						   -file   => $inputFile);
print "... done\n";

print "\nanalysing XML BLAST results\n";

# filtering threshold
my $threshold_coverage = 70;
my $threshold_expect = 0.00001;

$threshold_coverage = $input_cov if (defined $input_cov);
$threshold_expect = $input_exp if (defined $input_exp);

print "... coverage threshold: $threshold_coverage\n";
print "... expect threshold: $threshold_expect\n";

# remove .xml extension (if any)
$inputFile =~ s/\.xml$//;

my $o_file = "filter_".$inputFile;
open (OUT, ">$o_file") || die "cannot create output file!\n";

print "... exporting to file: $o_file\n";

print OUT "### $soft version $VERSION ($year)\n";
print OUT "### parameters: $param\n";
print OUT "### ".scalar localtime(time())."\n";
print OUT "#query\tq_length\thit\th_length\tHSP_length\tHSP_score(bits)\tHSP_expect\t%id\tq_coverage\tHit_coverage\n";

while( my $result = $in->next_result ) {
	#my $somethingToPrint = undef;
	my $q_name   = $result->query_description;
	my $q_length = $result->query_length;
	my $numHits  = $result->num_hits;
	
	
	while( my $hit = $result->next_hit ) {
		my $h_name   = $hit->name;
		my $h_length = $hit->length;

		
		while( my $hsp = $hit->next_hsp ) {
			my $hsp_length = $hsp->length('total');
			my $hsp_hitLength   = $hsp->length('hit');
			my $hsp_QueryLength = $hsp->length('query');
			my $hsp_expect = $hsp->evalue;
			my $hsp_bits   = $hsp->bits;
			my $percentID  = $hsp->percent_identity;
			
			my $coverage = ($hsp_QueryLength / $q_length) * 100;
			my $hit_coverage = ($hsp_hitLength / $h_length) * 100;
			
			next if ($hsp_expect > $threshold_expect);
			next if ($coverage < $threshold_coverage);
			
			print OUT "$q_name\t$q_length\t$h_name\t$h_length\t$hsp_length\t$hsp_bits\t$hsp_expect\t$percentID\t$coverage\t$hit_coverage\n";
			
			if (($hit_coverage > 100)||($coverage > 100)) {
				print "QUERY: $q_name - $q_length\nSUBJECT: $h_name - $h_length\n";
				print "HSP length: $hsp_length | Hit HSP length: $hsp_hitLength | Query HSP Length: $hsp_QueryLength\n";
				print "%id: $percentID | Query cov: $coverage\n";
			}
			

			
			#$somethingToPrint = defined;
		}
		
		
	}
	#print "--\n" if (defined $somethingToPrint);
}
print "... done\n";
close (OUT);

my $end_time=times();


print "\n--\n";
printf "Execution time: %.2f seconds CPU user-time\n",($end_time-$begin_time);