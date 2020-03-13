#!/usr/bin/perl
##This perl program captures and counts the adaptor ligated sequences at every protospacer position (Spacer Walking) and also the PAM 
##regions between the 5'flank and 3'flank (PAM Analysis)

##Usage: perl spacerWalkingAndPamAnalysis Experiment.txt 

##In the Spacer Walking mode, the Experiment.txt file(tab delimited) will contain all the protospacer positions and the associated flanks
##Example
## LibraryName	Position	Barcode	5_Flank	3_Flank
## Cas9	1	CTAGGT	CGGCATTCCTGCTGAACCGCTCTTCCGATCTA	AGTTGACCCA
## Cas9	2	CTAGGT	CGGCATTCCTGCTGAACCGCTCTTCCGATCTCA	AGTTGACCCA
## Cas9	3	CTAGGT	CGGCATTCCTGCTGAACCGCTCTTCCGATCTACA	AGTTGACCCA

##In the PAM Analysis mode, the Experiment.txt file will contain just the 1 position (highest activity) and its associated flanks 
##Example
## LibraryName	Position	Barcode	5_Flank	3_Flank
## Cas9	1	CTAGGT	CGGCATTCCTGCTGAACCGCTCTTCCGATCTACA	AGTTGACCCA

----------------------------------------------------------------------------------------------------------

##Define parameter values
my $inputFile = $ARGV[0];
my $barcodeLen=6;
my $pamLen=7;
my $bounds="both_flanks";
my $exp = "AK4130001";
my %flank_h={};
my $tsvdir = "TSV/";
mkdir $tsvdir;

open(F,$inputFile) || die "cannot open file"; 
while(<F>){
chomp;
my @arr = split("\t",$_);
$flank_h{$arr[2]}{$arr[1]}{"name"}=$arr[0];
$flank_h{$arr[2]}{$arr[1]}{"5_flank"}=$arr[3];
$flank_h{$arr[2]}{$arr[1]}{"3_flank"}=$arr[4];
}
close(F);
my @keys = keys %flank_h;
my @pos = sort keys %{ $flank_h{$keys[0]} };

##Loop over each position and count the number of PAM sequences between the 5_flank and 3_flank from the fastq file 
 foreach my $p(@pos){
	my $fname = $tsvdir.$exp."_".$p."_pam_expt.tsv"; 
	my $sname = $tsvdir.$exp."_".$p.".summary.txt"; 
	my $header = $exp."-".$p;

	my $tr=0; my $bct=0; my $mct=0; my %skpH=(); my %ctH=();

	open(FILE,">$fname");

	open(IN,"$exp.fastq") || die "cannot open file";
	$GET=0;
	while(<IN>){
		if(/^\@/o){$GET=1;next;} 
		if($GET){
			$tr++;$bc=substr($_,0,$barcodeLen);
			if( exists( $flank_h{$bc} ) ){
				if ($pamLen == 0) { exit; }
				$up=$flank_h{$bc}{$p}{"5_flank"};
				$dn=$flank_h{$bc}{$p}{"3_flank"};
				$MATCH = undef;
				if( $bounds eq "both_flanks" ) { if(/$up([ACGT]{$pamLen})$dn/){ $MATCH = $1;}} 
				elsif( $bounds eq "dn_flank" ) { if(/([ACGT]{$pamLen})$dn/){ $MATCH = $1;}}
				elsif( $bounds eq "up_flank" ) { if(/$up([ACGT]{$pamLen})/){ $MATCH = $1;}}
				elsif( $bounds eq "variable" ) { if(/$up([ACGT]+)$dn/){ $MATCH = $1; }}
				if( defined( $MATCH ) ) {
					$ctH{$bc}{$MATCH}+=1;
					$mct++;
				} else{chomp;$skpH{$bc}{$_}+=1;}
				$ctH{$bc}{"Total"}+=1;$bct++;}
		}
		$GET=0;
	}
	## Write the count file 
	foreach my $bc (sort (keys %ctH) ) {
		print FILE "#Expt\tSite\tCount\n";
		foreach $k (sort { $ctH{$bc}{$b}<=>$ctH{$bc}{$a} } (keys %{$ctH{$bc}}) ) {
			print FILE "$flank_h{$bc}{$p}{\"name\"}\t$k\t$ctH{$bc}{$k}\n";
		}
		print FILE "\nTop 50 rejected reads for $flank_h{$bc}{$p}{\"name\"}\n";
		$c=0;
		foreach $r (sort {$skpH{$bc}{$b}<=>$skpH{$bc}{$a}} (keys %{$skpH{$bc}})){
			print FILE "$r\t$skpH{$bc}{$r}\n";
			if($c>48){print FILE "\n";last;}
				$c++;
		}
	}

	close(FILE);
	
	## Write the summary file 
	open(S, ">$sname");
	print S "# $header $bounds\n# Total Reads:\t\t$tr\n# Barcode Matched:\t$bct\n# PAM Counted:\t\t$mct\n"; 
	close(S);

}
