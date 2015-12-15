#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# Name: lnc_02_analyzeBlastResults.pl
#-------------------------------------------------------------------------------
# Desc: Using filtered results from lnc_01_parseBlastResult, analyze TBLASTX
# results of predicted new GGA mRNA of FEELlnc vs HSA mRNA
#
# VERSION HISTORY:
# 1.0 : initial release
#-------------------------------------------------------------------------------
# Author: Frédéric lecerf on dec. 2015
#-------------------------------------------------------------------------------
use strict;
use Getopt::Long;

sub identity_rost ( $ );

#-------------------------------------------------------------------------------
# CONSTANT value definition
#-------------------------------------------------------------------------------
$|=1;
my $soft = "lnc_02_analyzeBlastResults";
my $VERSION = "1.0";
my $year = "2015";

# Benchmark 
my $begin_time = times();
my $verbose = undef;
my $blastFile = undef;

my @Options = @ARGV;

my %Options = (
        'blast=s'	 	=> \$blastFile,
		'verbose'		=> \$verbose,
	      );

my %OptionsHelp = (
		    'blast    ' => '-b [file], filtered XML blast output',
			'verbose  ' => '-v, verbose mode (optional)',
		    );

my $param = join(' ',@Options);
#-------------------------------------------------------------------------------
sub usage ( $ ) {
	my ($msg)=@_;
	print STDERR "Error: $msg\n";
	print "--------------------------------------------------------------------------------\n";
	print "$soft $VERSION ($year)\n";
	print "--------------------------------------------------------------------------------\n";
	print "Desc: Using filtered results from lnc_01_parseBlastResult, analyze TBLASTX\n";
	print "results of predicted new GGA mRNA of FEELlnc vs HSA mRNA\n";
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
defined $blastFile	|| &usage('filtered XML blast output required (-b)');


print "--------------------------------------------------------------------------------\n";
print "$soft $VERSION ($year)\n";
print "--------------------------------------------------------------------------------\n";

#-------------------------------------------------------------------------------
# Load blast filtered results
#-------------------------------------------------------------------------------

my %BlastResult = ();

my $nb_results = 0;

open (IN, $blastFile) || die "cannot open file: $blastFile\n";

print "Loading blast filered results from: $blastFile\n";
while (my $line = <IN>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#/);
	my @T = split (/\t/, $line);
	
	my $query      = $T[0];
	my $q_length   = $T[1];
	my $hit        = $T[2];
	my $h_length   = $T[3];
	my $HSP_length = $T[4];
	my $HSP_score  = $T[5];
	my $HSP_expect = $T[6];
	my $percentID  = $T[7];
	my $q_coverage = $T[8];
	my $h_coverage = $T[9];
	
	my @Query = split (/_/, $query);
	my @Hit   = split (/\|/, $hit);
	
	my $tcons = "TCONS".$Query[1];
	my $exons = $Query[6];
	my $enst = $Hit[0];
	my $ensg = $Hit[1];
	
	# compute various threshold (Rost & Li)
	# As TCONS are blasted by exons, the use of min is irrelevant here
	# We used only the R1 ratio to modify the %ID
	my $rost_threshold = identity_rost($HSP_length);
	my $r1 = $HSP_length / $q_length;
	#my $r2 = $HSP_length / $h_length;
	#my $min_ratio = $r1;
	#$min_ratio = $r2 if ($r2 < $r1);
	#my $I_prime = $percentID * $min_ratio;
	my $I_prime = $percentID * $r1;
	
	# Filtering results with the following conditions (3rd case in Li paper):
	my $I = 30;
	my $L = 150;
	# Apply filter
	if ($HSP_length >= $L) {
		next if $I_prime < $I;
	}
	else {
		next if $I_prime < $rost_threshold;
	}	
	
	$BlastResult{$tcons}{$exons}{$ensg}{$enst}{id} = $percentID;
	$BlastResult{$tcons}{$exons}{$ensg}{$enst}{idprime} = $I_prime;
	$BlastResult{$tcons}{$exons}{$ensg}{$enst}{qcov} = $q_coverage;
	$BlastResult{$tcons}{$exons}{$ensg}{$enst}{scov} = $h_coverage;

	$nb_results++;
	
	print "$tcons - $exons | $ensg\n";
}
close (IN);
print "... $nb_results Blast results loaded\n";
print "... ".(keys %BlastResult)." TCONS Transcripts loaded\n";
print "... done\n";





my $end_time=times();


print "\n--\n";
printf "Execution time: %.2f seconds CPU user-time\n",($end_time-$begin_time);



#-------------------------------------------------------------------------------
# Return ID threshold value depending on lenght (Rost, 1999)
#-------------------------------------------------------------------------------
sub identity_rost ( $ )
{
	my ($L) = @_;
	my $value = 0;
	if ($L >= 150)
	{
		$value = 30;
		return ($value);
	}
	my $n = 6;
	my $exposant = -0.32 * ( 1 + exp ( -($L/1000) ) );
	$value = ( 0.01*$n + 4.8 * ($L**$exposant) ) * 100;
	return ($value);
}