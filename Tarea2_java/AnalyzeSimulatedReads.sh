minOverlap=10

for i in SimulatedReads/*; 
    do
    echo ${i:15}
    java -Xmx4g -cp lib/NGSEPcore_3.2.0.jar:bin uniandes.algorithms.readsanalyzer.ReadsAnalyzerExample Kmers SimulatedReads/${i:15} 10 > temp
    timetable=$(grep "table" temp)
    echo -n $timetable" | "
    totalnumber=$(grep "Total number" temp)
    echo -n $totalnumber"  |"
    awk '/distribution/,0' temp | tail -n +2 > ProcsData/${i:15}"kmerSize_"10.xvg
    echo -n "mean abudance = "
    awk '{ sum += $2 } END { print(sum / NR) }' ProcsData/${i:15}"kmerSize_"10.xvg


    ### produce overlap-layout.

    java -Xmx4g -cp lib/NGSEPcore_3.2.0.jar:bin uniandes.algorithms.readsanalyzer.ReadsAnalyzerExample Overlap SimulatedReads/${i:15} $minOverlap > temp
    timetable=$(grep "overlap graph" temp) # print time
    echo -n $timetable" | "
    totalnumber=$(grep "Number of" temp) # print total number
    echo -n $totalnumber"  |"
    MeanAbud=$(grep "mean abundace" temp) # print mean abundance
    echo -n "$MeanAbud | "
    MeanSuc=$(grep "mean succesors" temp) # print mean succes
    echo -n "$MeanSuc | "
    MaxSuc=$(grep Max temp | tr -dc '0-9') # extract max number of succesors
    grep succ  -A$((MaxSuc +1)) temp | grep -v "succ" > ProcsData/${i:15}"Overlap_"$minOverlap.xvg # extract overlap distriburion in a file
    echo ""

    done