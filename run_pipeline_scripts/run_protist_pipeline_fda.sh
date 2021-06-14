#!/bin/sh
#SBATCH --job-name=taxaTarget
#SBATCH -t 7-0:00              # time limit: (D-HH:MM)
#SBATCH --mem=120G            # memory per node in MB
#SBATCH --nodes=1              # number of nodes
#SBATCH --cpus-per-task=12

# Load modules and software paths into environment
#
module load diamond
module load kaiju
module load biopython
module load python/3.8.1


# Set paths to directories and databases
#
createDB="/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/createDB_scripts/"
protist_data="/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/data/"
run_pipeline="/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/run_pipeline_scripts/"
kaijuDB=$protist_data"/marker_geneDB.fasta.kaiju.fmi"
queryDB=$protist_data"/marker_geneDB.fasta"


################################################################
################################################################
################################################################


echo taxaTarget $1

# Fastq file(s) to be analyzed
#
reads_fastq=$1

# Output file
#
out=$1.taxaTarget
mkdir $out
cd $out

# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
kaijux -f $kaijuDB -i ../$reads_fastq -z 12 -m 9 | grep "^C" > kaiju


# Extract reads that aligned to binning database
#
python $run_pipeline/extract_kaiju_reads.py -k kaiju -s ../$reads_fastq -o kaiju.fasta


# Align binned reads, with Diamond, to queryDB
#
# --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
diamond blastx --top 0 --sensitive --min-score 55 --db $queryDB --query kaiju.fasta --threads 12 --outfmt 6 --out kaiju.fasta.diamond


# Classify reads
#
python $run_pipeline/classify_reads_strict.py -n 0.95 -d kaiju.fasta.diamond -m $protist_data/marker_gene_metadata.txt -c $protist_data/strict_classifiers_filtered.txt -t 0.4 -f $protist_data/fullnamelineage.dmp -mg $protist_data/255_MGs.txt

# python $run_pipeline/classify_reads_multilogit.py -d kaiju.fasta.diamond -m $protist_data/marker_gene_metadata.txt -t 0.2 -p $protist_data/logistic_classifiers.txt -f $protist_data/fullnamelineage.dmp > reads_classified.txt
