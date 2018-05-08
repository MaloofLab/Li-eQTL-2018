#!/bin/bash
#SBATCH --partition=gc # partition to submit to
#SBATCH --job-name=“SNP_calling” # Job name
#SBATCH --array=1-40 # array
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=20 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem=128000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=00-10:00:00  # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=SNP_calling_01_24_2017_%A_%a.out # STDOUT
#SBATCH --error=SNP_calling_01_24_2017_%A_%a.err # STDERR
#SBATCH --mail-user=you@gmail.co # does not work yet
#SBATCH --mail-type=ALL # does not work yet
# This will be run once for a single process

/bin/hostname
module load gatk/3.6

# GATK doesn't support multithread so I will run chromosome intervals as array
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

THREADS=20

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`

reference=Brassica_napus_v4.1.chromosomes.fa

sample1=Ae_realignedBam.bam
sample2=Ol_realignedBam.bam

echo $reference
echo $sample1
echo $sample2
echo $sample 

# SNP calling for Ae and Ol combined by chromosome
# task array for differnet intervals
gatk    -T HaplotypeCaller\
        -R ${reference}\
        -I ${sample1}\
        -I ${sample2}\
        -dontUseSoftClippedBases\
        -stand_call_conf 20.0\
        -stand_emit_conf 20.0\
        -L ${sample}\
        -o output/${sample}.vcf

echo "done"
