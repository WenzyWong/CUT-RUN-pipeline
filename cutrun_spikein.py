######################################################
# CUT&RUN analysis for experiment with spike-in
# 
# Yunzhe Wang, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-08-28
# 
######################################################
import glob
import os

sample_dirs = glob.glob("Rawdata-clones/Panc-1-CUTRUN-*")
SAMPLE = [os.path.basename(d).replace("Panc-1-CUTRUN-", "") for d in sample_dirs]
INDEX_BT2 = "/data/yzwang/reference/hg38/hg38"
SPIKEIN_BT2 = "/data/yzwang/cowork/zxliu/cutrun/00_spk_ecoli_bt2/ecoli"
THREADS = 64

rule all:
    input:
        expand("01_spk_clean/{smpl}/{smpl}_R1.fq.gz", smpl = SAMPLE),
        expand("01_spk_clean/{smpl}/{smpl}_R2.fq.gz", smpl = SAMPLE),
        expand("02_spk_bam/{smpl}/{smpl}_bt2_hg38_sort.bam", smpl = SAMPLE),
        expand("02_spk_bam/{smpl}/{smpl}_bt2_hg38_sort.bam.bai", smpl = SAMPLE),
        expand("02_spk_bam_spikein/{smpl}/{smpl}_bt2_spikein_sort.bam", smpl = SAMPLE),
        expand("03_spk_bw_normalized/{smpl}/{smpl}_normalized_bin500.bw", smpl = SAMPLE),
        expand("04_spk_peak/{smpl}/{smpl}_peaks.broadPeak", smpl = SAMPLE),
        expand("04_spk_peak/{smpl}/{smpl}_peaks_adj.broadPeak", smpl = SAMPLE)

rule trim:
    input:
        "Rawdata-clones/Panc-1-CUTRUN-{smpl}/Panc-1-CUTRUN-{smpl}_R1.fq.gz",
        "Rawdata-clones/Panc-1-CUTRUN-{smpl}/Panc-1-CUTRUN-{smpl}_R2.fq.gz"
    output:
        "01_spk_clean/{smpl}/{smpl}_R1.fq.gz",
        "01_spk_clean/{smpl}/{smpl}_R2.fq.gz"
    log:
        "logs/01_spk_clean/{smpl}_cutadapt.log"
    shell:
        "cutadapt -j {THREADS} \
        -a AATGATACGGCGACCACCGAGATCTACAC \
        -A CAAGCAGAAGACGGCATACGAGAT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1"

rule bt2_mapping:
    input:
        "01_spk_clean/{smpl}/{smpl}_R1.fq.gz",
        "01_spk_clean/{smpl}/{smpl}_R2.fq.gz"
    output:
        "02_spk_bam/{smpl}/{smpl}_bt2_hg38.sam"
    log:
        "logs/02_spk_bam/{smpl}_bt2_hg38.log"
    shell:
        "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant \
        -x {INDEX_BT2} -p {THREADS} -1 {input[0]} -2 {input[1]} -S {output} > {log} 2>&1"

rule bt2_mapping_spikein:
    input:
        "01_spk_clean/{smpl}/{smpl}_R1.fq.gz",
        "01_spk_clean/{smpl}/{smpl}_R2.fq.gz"
    output:
        "02_spk_bam_spikein/{smpl}/{smpl}_bt2_spikein.sam"
    log:
        "logs/02_spk_bam/{smpl}_bt2_spikein.log"
    shell:
        "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant \
        -x {SPIKEIN_BT2} -p {THREADS} -1 {input[0]} -2 {input[1]} -S {output} > {log} 2>&1"

rule bam_file_sort:
    input:
        "02_spk_bam/{smpl}/{smpl}_bt2_hg38.sam"
    output:
        "02_spk_bam/{smpl}/{smpl}_bt2_hg38_sort.bam",
        "02_spk_bam/{smpl}/{smpl}_bt2_hg38_sort.bam.bai"
    log:
        "logs/02_spk_bam/{smpl}_bt2_hg38_sort.log"
    shell:
        """
        samtools sort -O BAM -o {output[0]} -T {output[0]}.temp -@ {THREADS} {input} 2>> {log}
        samtools index -@ {THREADS} {output[0]} {output[1]} 2>> {log}
        """

rule bam_file_sort_spikein:
    input:
        "02_spk_bam_spikein/{smpl}/{smpl}_bt2_spikein.sam"
    output:
        "02_spk_bam_spikein/{smpl}/{smpl}_bt2_spikein_sort.bam"
    log:
        "logs/02_spk_bam_spikein/{smpl}_bt2_spikein_sort.log"
    shell:
        "samtools sort -O BAM -o {output} -T {output}.temp -@ {THREADS} {input} > {log} 2>&1"

rule spikein_normalize_bin500:
    input:
        primary_bam = "02_spk_bam/{smpl}/{smpl}_bt2_hg38_sort.bam",
        spikein_bam = "02_spk_bam_spikein/{smpl}/{smpl}_bt2_spikein_sort.bam"
    output:
        "03_spk_bw_normalized/{smpl}/{smpl}_normalized_bin500.bw"
    params:
        outdir = "03_spk_bw_normalized/{smpl}"
    log:
        "logs/03_spk_bw_normalized/{smpl}_normalized_bin500.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        spikein_reads=$(samtools view -c -F 4 {input.spikein_bam})

        # paired-end: frag = reads / 2
        # scaling factor = 10000 / fragments = 10000 / (reads/2) = 20000 / reads
        if [ "$spikein_reads" -gt 0 ]; then
                scaling_factor=$(awk -v r="$spikein_reads" 'BEGIN{{printf "%.8f", 20000/r}}')
        else
            scaling_factor=1
        fi

        echo "Spike-in reads: $spikein_reads" >> {log}
        echo "Scaling factor: $scaling_factor" >> {log}

        bamCoverage --bam {input.primary_bam} -o {output} \
            --scaleFactor $scaling_factor \
            --binSize 500 --extendReads --normalizeUsing RPGC \
            --effectiveGenomeSize 2913022398 > {log} 2>&1
        """

rule call_peak:
    input:
        "02_spk_bam/{smpl}/{smpl}_bt2_hg38_sort.bam"
    output:
        "04_spk_peak/{smpl}/{smpl}_peaks.broadPeak"
    params:
        "{smpl}",
        "04_spk_peak/{smpl}"
    log:
        "logs/04_spk_peak/{smpl}_peaks.log"
    shell:
        "macs3 callpeak -q 0.1 -t {input} -g hs -n {params[0]} --broad --broad-cutoff 0.1 --outdir {params[1]} > {log} 2>&1"

rule adjust_peak:
    input:
        peak = "04_spk_peak/{smpl}/{smpl}_peaks.broadPeak",
        spikein_bam = "02_spk_bam_spikein/{smpl}/{smpl}_bt2_spikein_sort.bam"
    output:
        "04_spk_peak/{smpl}/{smpl}_peaks_adj.broadPeak"
    shell:
        r"""
        spikein_reads=$(samtools view -c -F 4 {input.spikein_bam})

        if [ "$spikein_reads" -gt 0 ]; then
            awk -v reads="$spikein_reads" 'BEGIN{{OFS="\t"; cf=20000/reads}} {{if(NR>1) $5=$5*cf; print}}' {input.peak} > {output}
        else
            cp {input.peak} {output}
        fi
        """