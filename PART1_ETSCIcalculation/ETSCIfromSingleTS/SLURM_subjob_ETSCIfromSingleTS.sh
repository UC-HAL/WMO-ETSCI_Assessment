#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --time=1-0:00:00
#SBATCH --mem=5G
#SBATCH --job-name=ETSCIfromSingleTS
#SBATCH --account=def-stadnykt-ab
#SBATCH --mail-user=alida.thiombiano@ucalgary.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#------------------------------------------------------ # Debugging section for PBS
echo "Node  $PBS_NODEFILE :"
echo "---------------------"
cat $PBS_NODEFILE
echo "---------------------"
echo "Shell is $SHELL"
NUM_PROCS=`/bin/awk 'END {print NR}' $PBS_NODEFILE`
echo "Running on $NUM_PROS processors."
echo "which mpirun = `which mpirun`"
#-------------------------------------------------------

echo `hostname`
echo "Current working directory is `pwd`"
echo "Starting run at: `date`"

# Write your commands below

#load modules required to run Climpact R software
module load r proj/7.2.1 udunits netcdf gdal geos

#copy as shortcuts (symbolic links), the files from the Climpact repo to the current working folder 
# with this, no need to specifiy the directory of the folder created after installing the Climpact software
ln -s ~/climpact/* ./

#Run the script to be executed
Rscript climpact.batch.stations.r ./inputSingleTS_txt/ ./metadata_batchSingleTS.txt 1985 2014 32

# end of your commands

#delete symbolic links (not used anymore)
find ./ -maxdepth 1 -type l -delete


echo "Parallel job finished with exit code $? at: `date`"
echo "****************************************************"

