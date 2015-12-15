# lnc_Analyses
some scripts to analyze Blast results from FEELnc prediction

---
lnc_01_parseBlastResults.pl

Desc: filter XML blast output (blastn, blastp, whatever) and return a tab delimited file with usefull informations (query and hit coeverage, name, length, e-value, %id).
Options:
	verbose  : -v, verbose mode (optional)
	input    : -i [file], XML blast output
	expect   : -e [nb], add expect threshold (optional, default 10-5)
	coverage : -c [nb], query coverage threshold value (in %, optional, default: 70%)

---
lnc_02_analyzeBlastResults.pl
Desc: Using filtered results from lnc_01_parseBlastResult, analyze TBLASTX results of predicted new GGA mRNA of FEELlnc vs HSA mRNA.
Options:
	verbose  : -v, verbose mode (optional)
	blast    : -b [file], filtered XML blast output
