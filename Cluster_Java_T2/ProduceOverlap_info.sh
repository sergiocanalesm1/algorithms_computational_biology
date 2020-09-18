# Generate overlap-Layout-Consensus in files of data directory
#
#
minOverlap=20
echo "Archivo         |  Tiempo (ms)|  Ab. Media | Media Suce. | Tiempo Ensamblaje |"
for i in data/te*; 
    do
    echo -n ${i:5}" | "
    java -Xmx4g -cp lib/NGSEPcore_3.2.0.jar:bin uniandes.algorithms.readsanalyzer.ReadsAnalyzerExample Overlap data/${i:5} $minOverlap > temp
    timetable=$(grep "overlap graph" temp | awk '{print $5}') # print time
    echo -n $timetable" | "
    totalnumber=$(grep "Number of" temp | awk '{print $5}') # print total number
   # echo -n $totalnumber"  |"
    MeanAbud=$(grep "mean abundace" temp | awk '{print $4}' ) # print mean abundance
    echo -n "$MeanAbud | "
    MeanSuc=$(grep "mean succesors" temp | awk '{print $4}') # print mean succes
    echo -n "$MeanSuc | "
    MaxSuc=$(grep Max temp | tr -dc '0-9') # extract max number of succesors
    grep succ  -A$((MaxSuc +1)) temp | grep -v "succ" > ProcsData/${i:5}"Overlap_"$minOverlap.xvg # extract overlap distriburion in a file
    assembl=$(grep "assembling" temp | awk '{print $4}' )
    echo -n "$assembl"
    echo ""
    #echo "Assembly"
    sed -e '1,/Assembly/ d' temp > Assembly/${i:5}"Assembly_Overlap_"$minOverlap.xvg
    done 
#rm temp 
