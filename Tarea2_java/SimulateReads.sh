# Simulate reads of sequencing process base on COVID.fasta obtained from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta
# varies seqeunece length in reads and also changes amount of secuences (5x)
Prof=(5 10 20 50 100)
SeqLength=(50 100 200 500)
#total bp in fasta file taken from http://seqanswers.com/forums/showthread.php?t=36160
TotalLeght=$(grep -v ">" $1 | wc | awk '{print $3-$1}') # total lenght of the sequence
error=0

if [ -z "$1" ]
  then
    echo "No fasta field found"
    exit 1
fi


if [ -z "$2" ]
  then
    echo "No error Found"
    exit 1
fi


 for prof in "${Prof[@]}"; 
        do  
        for seqLen in "${SeqLength[@]}"; 
              do
              TotalReads=$(( $prof*( $TotalLeght / $SeqLength ) ))  # calculate number of reads based on total lenght and seqlenght
              echo "Generating Simulated Read "$prof"X and seqlenght "$seqLen 
              java -Xmx4g -cp lib/NGSEPcore_3.2.0.jar:bin uniandes.algorithms.readsanalyzer.SimpleReadsSimulator $1 $seqLen $TotalReads SimulatedReads/$1$prof"X"_seqLen$seqLen"Error"$error".fastq"  $error
              done
        done
