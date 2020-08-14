#!/usr/bin/perl
use File::Basename;

## This perl program captures and counts the adaptor ligated sequences at every protospacer position (Spacer Walking) and also the PAM 
## regions between the 5'flank and 3'flank (PAM Analysis)

## Usage: perl spacerWalkingAndPamAnalysis RESULTS_DIR INPUT_FILE.txt EXPT.fastq BARCODE_LEN PAM_LEN
## Example: perl spacerWalkingAndPamAnalysis /path/to/resultsdir input.txt AK101.fastq 6 7
## ---------------------------------------------------------------------------------------------------------------------------------------------------

## Input Parameter values

my $dir = $ARGV[0];
my $inputFile = $ARGV[1];
my $fastq =$ARGV[2];
my $name = basename($fastq);
my $exp = $name;
$exp =~ s/.fastq//g;
my $barcodeLen=$ARGV[3];
my $pamLen=$ARGV[4];

my %flank_h={};
my $tsvdir = "$dir/PAM_Data/";

mkdir $tsvdir;

## Read the input file and store values in a hash

open(F,$inputFile) || die "cannot open file"; 
while(<F>){
chomp;
my @arr = split("\t",$_);
unless (/^Name/){
$flank_h{$arr[2]}{$arr[1]}{"name"}=$arr[0];
$flank_h{$arr[2]}{$arr[1]}{"5_flank"}=$arr[3];
$flank_h{$arr[2]}{$arr[1]}{"3_flank"}=$arr[4];
}
}
close(F);

my @keys = keys %flank_h;
my @pos = sort keys %{ $flank_h{$keys[0]} };

## Loop over each position and count the number of PAM sequences between the 5_flank and 3_flank from the fastq file 

foreach my $p(@pos){
	my $fname = $tsvdir.$exp."_".$p."_pam_expt.tsv"; 
	my $sname = $tsvdir.$exp."_".$p.".summary.txt"; 
	my $header = $exp."-".$p;

	my $tr=0; my $bct=0; my $mct=0; my %skpH=(); my %ctH=();

	open(FILE,">$fname");

	open(IN,"$fastq") || die "cannot open file";
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
				if(/$up([ACGT]{$pamLen})$dn/){ $MATCH = $1;}
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
		print FILE "#Sample\tPAM\tCount\n";
		foreach $k (sort { $ctH{$bc}{$b}<=>$ctH{$bc}{$a} } (keys %{$ctH{$bc}}) ) {
			print FILE "$flank_h{$bc}{$p}{\"name\"}\t$k\t$ctH{$bc}{$k}\n";
		}
		print FILE "\n\n";
	}

	close(FILE);
	
	## Write the summary file 

	open(S, ">$sname");
	print S "# $header\n# Total Reads:\t\t$tr\n# Barcode Matched:\t$bct\n# PAM Counted:\t\t$mct\n"; 
	close(S);

}

