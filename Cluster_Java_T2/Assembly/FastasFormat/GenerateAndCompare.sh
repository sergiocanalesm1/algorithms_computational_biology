### generate fasta files with lines of 70 b per line
for i in ../COVID*; do fold -w70 $i > ${i:3:25}.fasta ; done
#compare
 for i in COVID*; do echo $i; diff -U 0 ReferenceCovid.fa $i | grep -v ^@ | tail -n +3 | wc -l ; done

