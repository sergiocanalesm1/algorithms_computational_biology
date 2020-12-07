#!/bin/bash

# ###### Zona de Parámetros de solicitud de recursos a SLURM ############################
#
#SBATCH --job-name=Graphs_20frame_closeness		#Nombre del job
#SBATCH -p short			#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem=16000		#Memoria en Mb por CPU, Default=2048
#SBATCH --time=32:00:00			#Tiempo máximo de corrida, Default=2 horas
##SBATCH -o Graps_20frame_closenes.o%j			#Nombre de archivo de salida
#
####################################################################################

################# Zona Carga de Módulos ############################################
module load  anaconda/python3.7

####################################################################################


#### Zona de Ejecución de código y comandos a ejecutar secuencialmente #############
echo "Corri el:"
date
python main3.py 

echo -e "Finalicé la ejecución del script \n"
date
####################################################################################

