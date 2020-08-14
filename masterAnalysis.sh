#!/bin/bash

export configFile=$1

while read line
do
	if [[ "$line" =~ ^BARCODE.* ]];
	then
		barcodeLen=${line/BARCODE_LEN=/}
	fi
		if [[ "$line" =~ ^EXPECTED_PAM.* ]];
	then
		pamLen=${line/EXPECTED_PAM_LEN=/}
	fi
	if [[ "$line" =~ ^INPUT.* ]];
	then
		inputFile=${line/INPUT_FILE=/}
	fi
	if [[ "$line" =~ ^PROJECT.* ]];
	then
		project=${line/PROJECT_NAME=/}
	fi
	if [[ "$line" =~ ^EXPT.* ]];
	then
		expt=${line/EXPT_FASTQ=/}
	fi
	if [[ "$line" =~ ^RESULTS.* ]];
	then
		dir=${line/RESULTS_DIR=/}
	fi

done < $configFile

if [ ! -d $dir/$project ]; then
  	mkdir -p $dir/$project
	cp $configFile $dir/$project/config.txt
	cp $inputFile $dir/$project/input.txt
	ln -s $expt $dir/$project/
fi 

`/usr/bin/perl spacerWalkingAndPamAnalysis.pl $dir/$project $inputFile $expt $barcodeLen $pamLen`
#echo "/usr/bin/perl spacerWalkingAndPamAnalysis.pl $dir/$project $inputFile $expt $barcodeLen $pamLen"


