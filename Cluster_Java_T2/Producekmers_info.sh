# Generate kmers table of files of data directory
# varies kmersize around all files
# generate data of kmers distribution.

kmersize=(5 10 15 20 25 50 75)
for i in data/te*; 
    do
    echo $i
    for kmer in "${kmersize[@]}"; 
        do  
        echo -n " KmerSize = $kmer |"
        java -Xmx4g -cp lib/NGSEPcore_3.2.0.jar:bin uniandes.algorithms.readsanalyzer.ReadsAnalyzerExample Kmers data/${i:5} $kmer > temp
        timetable=$(grep "table" temp)
        echo -n $timetable" | "
        totalnumber=$(grep "Total number" temp)
        echo -n $totalnumber"  | "
        awk '/distribution/,0' temp | tail -n +2 > ProcsData/${i:5}"kmerSize_"$kmer.xvg
        echo -n "mean abudance = "
        awk '{ sum += $2 } END { print(sum / NR) }' ProcsData/${i:5}"kmerSize_"$kmer.xvg

     

        done
    done 
rm temp 
