En la carpeta del proyecto se encuentran los directories para la ejecuccion de programa en Java (bin,data,lib,src)
Adicionalmente se encuentran los directorios:

-Assembly/: Contiene los ensamblajes para los textos de ejemplo ( 100 anhos de soledad) y COVID ensamble. 
con formato. [NombreArchivo][Profundidad][seqLectura][error]Assembply[Overlap o kmer]
    .FastaFormat/: Directorio que contiene los archivos de ensamble en formato fasta con longitud de 70 b por lectura.

-ProceData/: Contiene los archivos de distribuciones de kmers y de overlap. 
con formato. [NombreArchivo][Profundidad][seqLectura][error][Overlap o kmer].xvg

    La longitud del genoma se estimo como : 
    awk '{print $1*$2}' COVID.fasta100X_seqLen200Error0.fastqkmerSize_20.xvg > sumaCOVID100X
    awk -F',' '{sum+=$1;} END{print sum/100;}' sumaCOVID100X

-SimulatedReads : Contiene los archivos fastq que salen de la simulacion de lecutras aleatoria.
con formato. [NombreArchivo][Profundidad][seqLectura][error].fastaq

-ErrorAnalisis: Contiene los archivos para la comparacion de las tasas de error con archivos 20x a 50 b por lectura.
-graphs: Contiene los archivos de las distribuciones.

Los datos reportados en el reporte fueron obtenidos al ejecutar los programas de Java sucesivamente mediante scripts de bash.

Producekmers_infor.sh :Produce una tabla de kmers para todos los archivos en data/ :
Muestra en consola una tabla donde se varia el kmer size y reporta tiempo,abundacia y numero total 
General los archivos: 
    [NombreArchivo][Profundidad][seqLectura][error]kmersize_[kmersize]
Se ejecuta en una terminal de linux como ./Producekmers_info.sh
Si se cambia data funciona para demas archivos fastq
Se puede cambiar el arrglo kmersize para explorar otros valores de kmersize en la ejecucion. 


ProduceOver_infor.sh : Produce grafos de sobrelapes para los archivos en data/
Muestra en consola tiempos de ejecucion y medias de abundacia y sucesores. 
General los archivos: 
    [NombreArchivo][Profundidad][seqLectura][error]Overlap_[minOverlap]
Se ejecuta en una terminal de linux como ./ProduceOver_info.sh
Si se cambia data funciona para otros archivos fastq

SimulateReads.sh: Produce los archivos de la carpeta SimulatedReads. 
    con formato. [NombreArchivo][Profundidad][seqLectura][error].fastaq
Se pueden variar los arreglos de profundidad y longitud de secuencia dentro de script para generar otros datos. 
Se ejecuta como ./SimulatedReads.sh [fastafile]. En la ejecucion se utilizo el archivo COVID.fasta


./AnalyzesSimulatedReads.sh 
Genera las distribuciones de kmers con tamaño de 20 y Grafos de sobrelapes para los archivos generados por SimulatedReads.sh


