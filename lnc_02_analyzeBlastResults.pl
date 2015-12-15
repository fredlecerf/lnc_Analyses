#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# Name: acg_05_analyzeBlastResults.pl
#-------------------------------------------------------------------------------
# Desc: Using filtered results from acg_04_parseBlastResult and the list of
# orthologous genes from acg_02_get_orthologsGenes, check if positive blast
# results are against a HSA or MMU orthologous genes present in the list
# Note : HSA and MMU analyzes have to be conducted separately with the
# -species parameter
#
# VERSION HISTORY:
# 1.0 : initial release
#-------------------------------------------------------------------------------
# Author: Frédéric lecerf on sept. 2013
#-------------------------------------------------------------------------------
use strict;
use Getopt::Long;


#-------------------------------------------------------------------------------
# CONSTANT value definition
#-------------------------------------------------------------------------------
$|=1;
my $soft = "acg_05_analyzeBlastResults.pl";
my $VERSION = "1.0";
my $year = "2013";

# Benchmark 
my $begin_time = times();
my $verbose = undef;
my $blastFile = undef;
my $orthoFile = undef;
my $species = undef;

my @Options = @ARGV;

my %Options = (
        'blast=s'	 	=> \$blastFile,
		'orthologs=s'	=> \$orthoFile,
		'species=s'		=> \$species,
		'verbose'		=> \$verbose,
	      );

my %OptionsHelp = (
		    'blast    ' => '-b [file], filtered XML blast output',
			'orthologs' => '-o [file], orthologous gene file',
			'species  ' => -'s [hsa or mmu], define ortholog species',
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
	print "Desc: Using filtered results from acg_04_parseBlastResult and the list of\n";
	print "orthologous genes from acg_02_get_orthologsGenes, check if positive blast\n";
	print "results are against a HSA or MMU orthologous genes present in the list\n";
	print "see source code description for more details\n";
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
defined $orthoFile	|| &usage('orthologous gene file required (-o)');
defined $species	|| &usage('acronym species required (-s)');


print "--------------------------------------------------------------------------------\n";
print "$soft $VERSION ($year)\n";
print "--------------------------------------------------------------------------------\n";

# check if species is correct
if (($species ne 'hsa')&&($species ne 'mmu')) {
	print STDERR "Usage: $0 [options]\n";
	print STDERR "see README for specific configuration parameters\n";
	print STDERR "Options:\n";
	map {printf STDERR "\t$_: %s\n",$OptionsHelp{$_};} keys %OptionsHelp; 
	print STDERR "Please cite: Lecerf F., $year\n";	
	print "========================================================================\n";
	die "Error: species name '$species' is not a valid acronym name\n";
}
else {
	print "species: $species\n";
}




#-------------------------------------------------------------------------------
# Load blast filtered results
#-------------------------------------------------------------------------------
my %Cuff_Trans2Genes = ();
my %Ens_Trans2Genes = ();
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
	
	my @Query = split (/\|/, $query);
	my @Hit   = split (/\|/, $hit);
	my $cuffGene = $Query[0];
	my $cuffTran = $Query[1];
	my $hitGene = $Hit[0];
	my $hitTran = $Hit[1];
	
	# conserve association between transcript and genes
	$Cuff_Trans2Genes{$cuffTran} = $cuffGene;
	$Ens_Trans2Genes{$hitTran} = $hitGene;
	
	# store blast results using only transcript ID
	$BlastResult{$cuffTran}{$hitTran}{q_length}   = $q_length;
	$BlastResult{$cuffTran}{$hitTran}{h_length}   = $h_length;
	$BlastResult{$cuffTran}{$hitTran}{HSP_length} = $HSP_length;
	$BlastResult{$cuffTran}{$hitTran}{HSP_score}  = $HSP_score;
	$BlastResult{$cuffTran}{$hitTran}{HSP_expect} = $HSP_expect;
	$BlastResult{$cuffTran}{$hitTran}{percentID}  = $percentID;
	$BlastResult{$cuffTran}{$hitTran}{q_coverage} = $q_coverage;
	$BlastResult{$cuffTran}{$hitTran}{h_coverage} = $h_coverage;
	$nb_results++;
}
close (IN);
print "... $nb_results Blast results loaded\n";
print "... ".(keys %BlastResult)." cufflinks Transcripts loaded\n";
print "... done\n";



#-------------------------------------------------------------------------------
# Load orthologous gene datafile
#-------------------------------------------------------------------------------
my %Orthologs = ();

open (IN, $orthoFile) || die "cannot open file: $orthoFile\n";

print "\nLoading orthologous genes data from: $orthoFile\n";
while (my $line = <IN>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#/);
	my @T = split (/\t/, $line);

	my $cuffgene = $T[0];
	my $cufftran = $T[1];
	
	my @OrthoGenes = undef;
	@OrthoGenes = split (',', $T[2]) if ($species eq 'hsa');
	@OrthoGenes = split (',', $T[3]) if ($species eq 'mmu');
	
	$Orthologs{$cufftran} = \@OrthoGenes;
}
close (IN);
print "... ".(keys %Orthologs)." cufflinks Transcript loaded\n";
print "... done\n";



#-------------------------------------------------------------------------------
# Analyze blast results taking account of orthologous genes
#-------------------------------------------------------------------------------
my $nb_noOrthologs = 0;
my $nb_HitWithOrthologs = 0;

print "\nanalyze blast results taking account of orthologous genes\n";
print "#cuffGene\tcuffTran\tEnsGene\thitTran\tq_length\th_length\tHSP_length\tHSP_score\tHSP_expect\tpercentID\tq_coverage\th_coverage\n";
foreach my $cuffTran (keys %BlastResult) {
	if (exists $Orthologs{$cuffTran}) {
		foreach my $hitTran (keys %{$BlastResult{$cuffTran}}) {
			my $EnsGene  = $Ens_Trans2Genes{$hitTran};
			my $cuffGene = $Cuff_Trans2Genes{$cuffTran};
			
			if ($EnsGene ~~ @{$Orthologs{$cuffTran}}) {
				my $q_length = $BlastResult{$cuffTran}{$hitTran}{q_length};
				my $h_length = $BlastResult{$cuffTran}{$hitTran}{h_length};
				my $HSP_length = $BlastResult{$cuffTran}{$hitTran}{HSP_length};
				my $HSP_score  = $BlastResult{$cuffTran}{$hitTran}{HSP_score};
				my $HSP_expect = $BlastResult{$cuffTran}{$hitTran}{HSP_expect};
				my $percentID  = $BlastResult{$cuffTran}{$hitTran}{percentID};
				my $q_coverage = $BlastResult{$cuffTran}{$hitTran}{q_coverage};
				my $h_coverage = $BlastResult{$cuffTran}{$hitTran}{h_coverage};
				
				print "$cuffGene\t$cuffTran\t$EnsGene\t$hitTran\t$q_length\t$h_length\t$HSP_length\t$HSP_score\t$HSP_expect\t$percentID\t$q_coverage\t$h_coverage\n";
				$nb_HitWithOrthologs++;
			}
			
		}
	}
	else {
		$nb_noOrthologs++;
	}
	
	
}
print "... $nb_noOrthologs cufflinks transcripts with no orthologous genes in $species\n";
print "... $nb_HitWithOrthologs Blast hit of cufflinks transcripts with orthologs\n";
print "... done\n";



my $end_time=times();


print "\n--\n";
printf "Execution time: %.2f seconds CPU user-time\n",($end_time-$begin_time);