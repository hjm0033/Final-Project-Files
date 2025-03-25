#!/bin/sh
 
######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to 
##    Use HiSat2 to index your reference genome and then map your cleaned (paired) reads to the indexed reference. If you have a large genome, this will require a lot more resources then listed here, such as the large queue with 1 core and 120gb. You may want to do the indexing as a seperate script.
##              First need to use gffread to convert annotation file from .gff3 to .gft format.
##              Use Stringtie to count the reads mapped to genes and transcripts, defined in this case by the genome annotation file
##              use the python script to take the Stringtie results to make two counts matricies, one at the gene level and one at the transcript level
## HiSat2  Indexing   InPut: Reference genome file (.fasta), and annotation file (.gff3) (Optional)
##                    Output: Indexed genome 
## HiSat2 Mapping     Input: Cleaned read files, paired (.fasq); Indexed genome
##                    Output: Alignment .sam files  
## Samtools  Convert .sam to .bam and sort          Input: Alignment files,  .sam
##                                                  Output: Sorted  .bam files
## Stringtie  Counting reads  Input: sorted .bam file
##                            Output:  Directories of counts files for Ballgown (R program for DGE)
##              prepDE.py    Python script to create a counts matrics from the Stringtie output.  Inputs: Directory from Stringtie
##                                                                                                Output:  .csv files of counts matrix
## For running the script on the Alabama Super Computer.
##  For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
##  After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
##  then run the script by using "run_script [script name]"
##  suggested paramenters are below to submit this script.
##    queue: class or medium
##    core: 6
##    time limit (HH:MM:SS): 04:00:00 
##    Memory: 12gb
##    
###############################################

#### Load all the programs you are going to use in this script.
source /apps/profiles/modules_asax.sh.dyn
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools/1.18
module load bcftools/1.13
module load gffread/0.12.7
#module load gffcompare


#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
  ## Replace the [#] with paths to define these variable
MyID=aubclsd0316          ## Example: MyID=aubtss

WD=/scratch/$MyID/DogBrainProjectTemp                      ## Example:/scratch/$MyID/PracticeRNAseq  
CD=/scratch/$MyID/DogBrainProjectTemp/CleanData             ## Example:/scratch/$MyID/PracticeRNAseq/CleanData   #   *** This is where the cleaned paired files are located from the last script
REFD=/scratch/$MyID/DogBrainProjectTemp/TashaRefGenome    ## Example:/scratch/$MyID/PracticeRNAseq/DaphniaRefGenome    # this directory contains the indexed reference genome for the garter snake
MAPD=/scratch/$MyID/DogBrainProjectTemp/Map_HiSat2           ## Example:/scratch/$MyID/PracticeRNAseq/Map_HiSat2      #
COUNTSD=/scratch/$MyID/DogBrainProjectTemp/Counts_StringTie       ## Example:/scratch/$MyID/PracticeRNAseq/Counts_StringTie
RESULTSD=/home/$MyID/DogBrainProjectTemp/Counts_H_S          ## Example:/home/aubtss/PracticeRNAseq/Counts_H_S

REF=GCF_000002285.5_Dog10K_Boxer_Tasha_genomic                   ## This is what the "easy name" will be for the genome reference

## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD
### Copy the reference genome (.fasta) and the annotation file (.gff3) to this REFD directory
cp /home/${MyID}/DogBrainProject/ncbi_dataset/data/GCF_000002285.5/${REF}.fna .
cp /home/${MyID}/DogBrainProject/ncbi_dataset/data/GCF_000002285.5/genomic.gff .

###  Identify exons and splice sites on the reference genome
gffread genomic.gff -T -o ${REF}.gtf               ## gffread converts the annotation file from .gff3 to .gft formate for HiSat2 to use.
hisat2_extract_splice_sites.py ${REF}.gtf > ${REF}.ss
hisat2_extract_exons.py ${REF}.gtf > ${REF}.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
hisat2-build --ss ${REF}.ss --exon ${REF}.exon ${REF}.fna Tasha_index

