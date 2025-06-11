######################################################
# CUT&RUN analysis for pre-experiment with no spike-in
# 
# Yunzhe Wang, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-06-11
# 
######################################################

SAMPLE = {"H3K9me3-1", "IgG"}
INDEX_BT2 = "/data/yzwang/reference/hg38/hg38"
THREADS = 80

rule all:
    input:
        expand("01_clean/{smpl}/{smpl}_R1.fq.gz", smpl = SAMPLE),
        expand("01_clean/{smpl}/{smpl}_R2.fq.gz", smpl = SAMPLE),
        expand("02_bam/{smpl}/{smpl}_bt2_hg38_sort.bam", smpl = SAMPLE),
        expand("02_bam/{smpl}/{smpl}_bt2_hg38_sort.bam.bai", smpl = SAMPLE),
        expand("03_bw/{smpl}/{smpl}_bt2_hg38.bw", smpl = SAMPLE),
        expand("03_bw/{smpl}/{smpl}_bt2_hg38_rpgc_bin500.bw", smpl = SAMPLE),
        expand("04_peak/{smpl}/{smpl}_peaks.broadPeak", smpl = SAMPLE)

rule trim:
    input:
        "Rawdata/Panc-1-CUT-RUN-{smpl}/Panc-1-CUT-RUN-{smpl}_R1.fq.gz",
        "Rawdata/Panc-1-CUT-RUN-{smpl}/Panc-1-CUT-RUN-{smpl}_R2.fq.gz"
    output:
        "01_clean/{smpl}/{smpl}_R1.fq.gz",
        "01_clean/{smpl}/{smpl}_R2.fq.gz"
    log:
        "logs/01_clean/{smpl}_cutadapt.log"
    shell:
        "cutadapt -j {THREADS} \
        -a AATGATACGGCGACCACCGAGATCTACAC \
        -A CAAGCAGAAGACGGCATACGAGAT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1"

rule bt2_mapping:
    input:
        "01_clean/{smpl}/{smpl}_R1.fq.gz",
        "01_clean/{smpl}/{smpl}_R2.fq.gz"
    output:
        "02_bam/{smpl}/{smpl}_bt2_hg38.sam"
    log:
        "logs/02_bam/{smpl}_bt2_hg38.log"
    shell:
        "bowtie2  --end-to-end --very-sensitive --no-mixed --no-discordant \
        -x {INDEX_BT2} -p {THREADS} -1 {input[0]} -2 {input[1]} -S {output} > {log} 2>&1"

rule bam_file_sort:
    input:
        "02_bam/{smpl}/{smpl}_bt2_hg38.sam"
    output:
        "02_bam/{smpl}/{smpl}_bt2_hg38_sort.bam",
        "02_bam/{smpl}/{smpl}_bt2_hg38_sort.bam.bai"
    log:
        "logs/02_bam/{smpl}_bt2_hg38_sort.log"
    shell:
        """
        samtools sort -O BAM -o {output[0]} -T {output[0]}.temp -@ {THREADS} {input} 2>> {log}
        samtools index -@ {THREADS} {output[0]} {output[1]} 2>> {log}
        """

rule bam_to_bw_rpgc:
    input:
        "02_bam/{smpl}/{smpl}_bt2_hg38_sort.bam"
    output:
        "03_bw/{smpl}/{smpl}_bt2_hg38_rpgc_bin500.bw"
    log:
        "logs/03_bw/{smpl}_bw_rpgc_bin500.log"
    shell:
        "bamCoverage --bam {input} -o {output} -bs 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --binSize 500 --extendReads > {log} 2>&1"

rule call_peak:
    input:
        "02_bam/{smpl}/{smpl}_bt2_hg38_sort.bam"
    output:
        "04_peak/{smpl}/{smpl}_peaks.broadPeak"
    params:
        "{smpl}",
        "04_peak/{smpl}"
    log:
        "logs/04_peak/{smpl}_peaks.log"
    shell:
        "macs3 callpeak -q 0.1 -t {input} -g hs -n {params[0]} --broad --broad-cutoff 0.1 --outdir {params[1]} > {log} 2>&1"
