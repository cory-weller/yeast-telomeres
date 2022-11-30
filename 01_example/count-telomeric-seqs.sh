#!/usr/bin/env bash
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=5:00:00
#SBATCH --gres=lscratch:120
#SBATCH --partition=norm
#SBATCH --output %j.slurm.out
#SBATCH --error %j.slurm.out

SRA=${1}
MINREPEATS=(5 10 20 30 35)
DATADIR='/home/wellerca/pacbio-yeast-genomes/data/'

cd /lscratch/${SLURM_JOB_ID}

unzip "${DATADIR}/shortreads/${SRA}.fastq.zip"                     # Inflate all fastqs


wc -l *.fastq > ${SRA}_linecounts.txt                              # count lines

for i in ${MINREPEATS[@]}; do
    pattern="(CA|CCA|CCCA){${i},}|(TG|TGG|TGGG){${i},}"
    grep -c -E ${pattern} *.fastq > ${SRA}_${i}_telo-counts.txt    # count telomere seqs
done


mv *.txt ${DATADIR}

rm *.fastq 





# for file in *.fastq; do                                     # Iterate over fastqs
#     echo ${file}
#     nLines=$(wc -l ${file} | awk '{print $1}')                  # Count all lines in file
#     nReads=$(python -c "print(${nLines}/4)")                    # Report nReads as NLines/4
#     nTelos=$(grep -c -E ${pattern} ${file})                     # Count telomeric reads
#     echo -e ${file%.fastq}\\t${nReads}\\t${nTelos} >> ${SRA}_telo-counts.txt
# done

