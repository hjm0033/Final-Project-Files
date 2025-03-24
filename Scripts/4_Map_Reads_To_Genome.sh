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
module load samtools
module load bcftools
module load gffread
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

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

# Move to the data directory
cd ${CD}  #### This is where our clean paired reads are located.

## Create list of fastq files to map.    Example file format of your cleaned reads file names: SRR629651_1_paired.fastq SRR629651_2_paired.fastq
## grab all fastq files, cut on the underscore, use only the first of the cuts, sort, use unique put in list
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list2    #should list Example: SRR629651

## Move to the directory for mapping
cd ${MAPD}

## move the list of unique ids from the original files to map
mv ${CD}/list2 .

## Process the samples in the list one by one using a while loop
while read i;
do
  ## Define file paths for paired and unpaired reads
  PAIRED_1="${CD}/${i}_1_paired.fastq"
  PAIRED_2="${CD}/${i}_2_paired.fastq"
  UNPAIRED="${CD}/${i}_unpaired.fastq"

  ## HiSat2 is the mapping program
  ## -p indicates number of processors, --dta reports alignments for StringTie
  ## Check if paired-end files exist, otherwise run in single-end mode
  if [[ -f "$PAIRED_1" && -f "$PAIRED_2" ]]; then
    echo "Running HiSat2 in paired-end mode for sample: $i"
    hisat2 -p 6 --dta --phred33 \
      -x "${REFD}/Tasha_index" \
      -1 "$PAIRED_1" -2 "$PAIRED_2" \
      -S "${i}.sam"

  elif [[ -f "$UNPAIRED" ]]; then
    echo "Running HiSat2 in single-end mode for sample: $i"
    hisat2 -p 6 --dta --phred33 \
      -x "${REFD}/Tasha_index" \
      -U "$UNPAIRED" \
      -S "${i}.sam"

  else
    echo "Warning: No valid FASTQ files found for sample: $i"
  fi

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
  samtools view -@ 6 -bS "$i".sam > "$i".bam  

    ###  This is sorting the bam, using 6 threads, and producing a .bam file that includes the word 'sorted' in the name
  samtools sort -@ 6  "$i".bam  -o  "$i"_sorted.bam

    ### Index the BAM and get mapping statistics, and put them in a text file for us to look at.
  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

  ### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model. 
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to 
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 
  
mkdir "${COUNTSD}"/"$i"
stringtie -p 6 -e -B -G  "${REFD}"/"${REF}".gtf -o "${COUNTSD}"/"$i"/"$i".gtf -l "$i"   "${MAPD}"/"$i"_sorted.bam

done<list2

#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *.txt ${RESULTSD}

### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix. 
 ## Move to the counts directory
cd ${COUNTSD}
 ## run the python script prepDE.phy to prepare you data for downstream analysis.
cp /home/${MyID}/class_shared/prepDE.py3 .

 prepDE.py3 -i ${COUNTSD}

### copy the final results files (the count matricies that are .cvs) to your home directory. 
cp *.csv ${RESULTSD}
## move these results files to your personal computer for downstream statistical analyses in R studio.
