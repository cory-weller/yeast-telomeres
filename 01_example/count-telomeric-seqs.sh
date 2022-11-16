#!/usr/bin/env bash
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --gres=lscratch:120
#SBATCH --partition=quick
#SBATCH --output %j.slurm.out
#SBATCH --error %j.slurm.out

SRA=${1}

DATADIR='/home/wellerca/pacbio-yeast-genomes/data'
pattern="(CA|CCA|CCCA){8,}|(TG|TGG|TGGG){8,}"

cd /lscratch/${SLURM_JOB_ID}

echo -e "SRA\\tnReads\\tnTeloReads" > telo-counts.txt

unzip "${DATADIR}/shortreads/${SRA}.fastq.zip"              # Inflate all fastqs

for file in *.fastq; do                                     # Iterate over fastqs
    echo ${file}
    nLines=$(wc -l ${file} | awk '{print $1}')                  # Count all lines in file
    nReads=$(python -c "print(${nLines}/4)")                    # Report nReads as NLines/4
    nTelos=$(grep -c -E ${pattern} ${file})                     # Count telomeric reads
    echo -e ${file%.fastq}\\t${nReads}\\t${nTelos} >> ${SRA}_telo-counts.txt
done

mv ${SRA}_telo-counts.txt ${DATADIR}/${SRA}_telo-counts.txt

rm *.fastq 

cd