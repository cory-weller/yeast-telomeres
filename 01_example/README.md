# Yeast Telomere Read Counts

For a first pass, do crude check of 'telomere-like reads' for all individuals in crosses containing YJM981.

| Strain1 | Strain2 | SRA run ID |
|:--------|:--------|:-----------|
| YJM981  | CBS2888 | SRR9330808 |
| CLIB219 | CBS2888 | SRR9330831 |
| YJM981  | 273614  | SRR9330809 |
| PW5     | 273614  | SRR9330830 |

[`count-telomeric-seqs.sh`](count-telomeric-seqs.sh) takes a single argument, the SRA run of interest, and 
unzips all relevent `fastq` files to an `lscratch` disk allocation. Then, it counts the number of reads
within a sample's `fastq` file that match a regular expression crudely identifying telomere-like sequences.
The regular expression is

```bash
"(CA|CCA|CCCA){8,}|(TG|TGG|TGGG){8,}"
```

which matches 8+ repeats of the telomeric unit `[C,CC,CCC]A` or `T[G,GG,GGG]` then taking the fraction of reads that are telomere-like multiplied by 1M to get RPM (reads per million).

Plotting the resulting data with [`plot-RPM-distributions.R`](plot-RPM-distributions.R)

![](telomere-rpm-01.png)

The crosses with YJM981 seem to be bimodal, but the additional peak appears to be higher than for crosses without YJM981.

Then let's just go ahead and count for all crosses:
```bash
sbatch count-telomeric-seqs.sh SRR9330808
sbatch count-telomeric-seqs.sh SRR9330831
sbatch count-telomeric-seqs.sh SRR9330809
sbatch count-telomeric-seqs.sh SRR9330830
sbatch count-telomeric-seqs.sh SRR9330832
sbatch count-telomeric-seqs.sh SRR9330811
sbatch count-telomeric-seqs.sh SRR9330837
sbatch count-telomeric-seqs.sh SRR9330810
sbatch count-telomeric-seqs.sh SRR9330836
sbatch count-telomeric-seqs.sh SRR9330813
sbatch count-telomeric-seqs.sh SRR9330815
sbatch count-telomeric-seqs.sh SRR9330812
sbatch count-telomeric-seqs.sh SRR9330833
sbatch count-telomeric-seqs.sh SRR9330814
sbatch count-telomeric-seqs.sh SRR9330817
sbatch count-telomeric-seqs.sh SRR9330816
```

![](telomeric-min8x.png)

${STRAIN}_${CHR}.fasta

sed "s/${match}/LOL/g" <( cat /data/SBGE/cory/pacbio/assemblies/CBS2888_chrI.fasta | tr -d '\n') | fold -w 80 | head

initial stretch of AC nucleotides, and final stretch of GT nucleotides

strain="273614"
chromosomes="chr{I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI,}"

match=$(cat /data/SBGE/cory/pacbio/assemblies/CBS2888_chrI.fasta | \
    awk 'NR > 1' | \
    head -n 500 | \
    tr -d '\n' | \
    grep -o -E "^[AC]{2,}")

cat /data/SBGE/cory/pacbio/assemblies/CBS2888_chrI.fasta | \
    tail -n 500 | \
    tr -d '\n' | \
    grep -o -E "[GT]{2,}$"


## Split up genomic and telomeric sequences

Using the `remove_telomeres` function, for each strain I generated two `fasta` files: one that
contains the leading/terminal telomeric sequences, and one that contains everything between those
telomeric sequences.

```bash
remove_telomeres() {
    local filename=${1}
    local strain=$(echo ${filename} | cut -d '_' -f 1)
    local chr=$(echo ${filename%.fasta} | cut -d '_' -f 2)
    echo ${filename}

    local L_pattern="^.*?(?=[TG][TG]|[TG].[TG])"
    local L_telomere=$(cat ${filename}    | \
                    awk 'NR > 1' | \
                    head -n 20 | \
                    tr -d '\n' | \
                    grep -o -P ${L_pattern} )

    local R_pattern="^.*?(?=[AC][AC]|[AC].[AC])"
    local R_telomere=$(cat ${filename}    | \
                    awk 'NR > 1' | \
                tail -n 20    | \
                tr -d '\n'      | \
                rev | \
                grep -o -P ${R_pattern} | \
                rev)
    
    # write telomeres to combined (strain-specific) fasta file
    echo ">${strain}_${chr}_L" >> ${strain}_telomeres.fasta
    echo ${L_telomere} | fold -w 60 >> ${strain}_telomeres.fasta
    echo ">${strain}_${chr}_R" >> ${strain}_telomeres.fasta
    echo ${R_telomere} | fold -w 60 >> ${strain}_telomeres.fasta

    # generate combined (strain-sepcific) fasta excluding terminal telomeric sequences
    echo ">${strain}_${chr}_minus_telomere" >> ${strain}_minus_telomeres.fasta
    awk 'NR > 1' ${filename}        | \
        tr -d '\n'                  | \
        sed "s/^${L_telomere}//g"   | \
        sed "s/${R_telomere}$//g"   | \
        fold -w 60 >> ${strain}_minus_telomeres.fasta
}

# While in directory containing ${STRAIN}_${CHR}.fasta files:
for file in *.fasta; do
    remove_telomeres ${file}
done

# Generate zip file for newly made fastas
zip 16strains_telomere_fasta *telomere*
```

## Soak up non-telomeric reads

The `${STRAIN}_minus_telomeres.fasta` file serves as the reference for `bwa` read mapping.

```bash
get_unmapped_reads() {
    local strain1=${1}
    local strain2=${2}
    local srr=${3}

    # Modules
    module load bwa/0.7.17
    module load samtools/1.16.1

    # move to temporary directory
    local datadir="/data/SBGE/cory/pacbio/assemblies"
    local gitdir="/home/wellerca/yeast-telomeres/01_example"

    local tmpdir="/lscratch/${SLURM_JOB_ID}"
    mkdir -p ${tmpdir} && cd ${tmpdir}


    # unpack fastqs
    local fastqs=${datadir}/${shortreads}/${srr}.fastq.zip
    unzip ${fastqs}

    # unpack ref genomes and concatenate into single reference
    # to soak up reads mapping to either parent in the cross
    unzip -p ${gitdir}/16strains_telomere_fasta.zip ${strain1}_minus_telomeres.fasta >  ref.fasta
    unzip -p ${gitdir}/16strains_telomere_fasta.zip ${strain2}_minus_telomeres.fasta >> ref.fasta

    # index concatenated ref
    bwa index ref.fasta

    # READY FOR SOAKING
    for file in *.fastq; do
        local sample="${file%.fastq}"

        bwa mem ${file} ref.fasta | \
            samtools view -hb - | \
            samtools sort -f 4 - > ${sample}_unmapped.bam
        samtools fastq -o ${sample}_unmapped.fastq ${sample}_unmapped.bam
    done
}


273614 YJM981 SRR9330809 3003
YJM981 CBS2888  SRR9330808 3004


    # alternative to samtools fastq
    # module load bedtools
    # bedtools bamtofastq -i ${sample}.bam -fq ${srr}_${sample}.fastq

    # extract unmapped reads

}

get_unmapped_reads YJM978_minus_telomeres.fasta /data/SBGE/cory/pacbio/shortreads/YJM978
```

map_and_call() {
    local zipfile=${1}
    local fastqfile=${2}
    local id=${fastqfile%.fastq}

    bwa mem mask.fasta <((unzip -p ${zipfile} ${fastqfile})) | \
        samtools view -hb - | \
        samtools sort - > ${id}.bam

    samtools index ${id}.bam
    bcftools mpileup \
        --targets-file targets.txt \
        --fasta-ref mask.fasta \
        ${id}.bam | \
        bcftools call \
        --ploidy 1 -m -Ob | \
        bcftools view | \
        sed 's/1:.*$/1/g' | \
        grep -v "^##" | awk '{print $1,$2,$4,$5,$10}' | sed 's/\.bam$//g' > ${id}.call
    
    rm ${id}.{bam,bam.bai}
}
export -f map_and_call



```bash
# Returns starts and ends of fasta file, up until the seq no longer looks 'telomeric'
# e.g. until 5' end has 2 [TG] nucleotides within a triple
# or   until 3' end has 2 [AC] nucleotides within a triple
test_telos() {
    local filename=${1}
    local L_pattern="^.*?(?=[TG][TG]|[TG].[TG])"
    local R_pattern="^.*?(?=[AC][AC]|[AC].[AC])"
    local L_telomere=$(cat ${filename}    | \
                    awk 'NR > 1' | \
                    head -n 20 | \
                    tr -d '\n' | \
                    grep -o -P ${L_pattern} )
    local R_telomere=$(cat ${filename}    | \
                    awk 'NR > 1' | \
                tail -n 20    | \
                tr -d '\n'      | \
                rev | \
                grep -o -P ${R_pattern} | \
                rev)

    echo ${L_telomere}
    echo ""
    echo ${R_telomere}
}

test_telos test.fasta
test_telos bad.fasta
```