# taxaTarget
taxaTarget is a tool for the classification of eukaryotes from metagenomic reads. The input is a fastq file of raw metagenomic sequencing reads. The reads are first mapped to the database with Kaiju to quickly identify reads that are likely sequenced from eukaryotic marker genes. The binned reads are then more sensitively aligned with Diamond to the marker genes. The Diamond output is then provided as input to classifiers, which were trained using the UniRef100 database and a supervised learning approach. Finally, the results are aggregated together and the taxonomic profile of the metagenomic sample is generated.

# Requirements
* Bullet list
Python3 (version 3.6 or higher; was tested on 3.7)
Pandas
Numpy
 
