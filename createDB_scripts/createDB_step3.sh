#!/bin/sh
#SBATCH --job-name=busco
#SBATCH -t 7-0:00              # time limit: (D-HH:MM)
#SBATCH --mem=120G            # memory per node in MB
#SBATCH --nodes=1              # number of nodes
#SBATCH --cpus-per-task=12


# Load modules and software paths into environment
#
module load biopython
module load diamond
module load kaiju


PWD=`pwd`
createDB=$PWD
data=$PWD"/../data/"
busco=$PWD"/../busco/"
run_pipeline=$PWD"/../run_pipeline_scripts/"


# Change to data directory
#
mkdir $data
cd $data


# Retrieve classification thresholds for each protein
#
python $createDB/strict_threshold_classifiers.py -d $data/diamond_results/ -m $data/marker_gene_metadata.txt


# Collect summary information about the classification power of marker genes
#
python $createDB/classification_power_per_protein_region.py



# Filter regions for species level classification power
#

