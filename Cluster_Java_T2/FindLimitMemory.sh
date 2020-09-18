#!/bin/bash

f=COVID.fasta
#if [ -z "$1" ]; then
#    echo " Missing fasta file"
#    exit 1
#fi
Prof=(120000 130000 110000 140000 )
TotalLeght=$(grep -v ">" $f | wc | awk '{print $3-$1}') # total lenght of the sequence
error=0
seqLen=50

minOverlap=20
for prof in "${Prof[@]}";
        do 
        echo "Evaluating $prof" 
 
        TotalReads=$(($prof*( $TotalLeght / $seqLen )))  # calculate number of reads based on total lenght and seqlenght
	java -Xmx4g -cp lib/NGSEPcore_3.2.0.jar:bin uniandes.algorithms.readsanalyzer.SimpleReadsSimulator $f $seqLen $TotalReads LimitMemory/$f$prof"X"_seqLen$seqLen"Error"$error".fastq"  $error
	java -Xmx4g -cp lib/NGSEPcore_3.2.0.jar:bin uniandes.algorithms.readsanalyzer.ReadsAnalyzerExample Overlap LimitMemory/$f$prof"X"_seqLen$seqLen"Error"$error".fastq" $minOverlap > temp
        timetable=$(grep "overlap graph" temp) # print time
        echo -n $timetable" | "
        totalnumber=$(grep "Number of" temp) # print total number
        echo -n $totalnumber"  |"
        MeanAbud=$(grep "mean abundace" temp) # print mean abundance
        echo -n "$MeanAbud | "
        MeanSuc=$(grep "mean succesors" temp) # print mean succes
        echo -n "$MeanSuc | "
        MaxSuc=$(grep Max temp | tr -dc '0-9') # extract max number of succesors
        grep succ  -A$((MaxSuc +1)) temp | grep -v "succ" > LimitMemory/${i:15}"Overlap_"$minOverlap.xvg # extract overlap distriburion in a file
        echo ""
        sed -e '1,/Assembly/ d' temp > LimitMemory/$f$prof"X"_seqLen$seqLen"Error"$error"_Assembly_Overlap_"$minOverlap.xvg


done
