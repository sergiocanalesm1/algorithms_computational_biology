### generate fasta files with lines of 70 b per line
for i in ../COVID*; do fold -w70 $i > ${i:3:25}.fasta ; done


