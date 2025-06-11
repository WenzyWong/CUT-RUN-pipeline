sortbam_k9="02_bam/H3K9me3-1/H3K9me3-1_bt2_hg38_sort.bam"
sortbam_igg="02_bam/IgG/IgG_bt2_hg38_sort.bam"

bw_k9="03_bw/H3K9me3-1/H3K9me3-1_bt2_hg38_rpgc_bin500.bw"
bw_igg="03_bw/IgG/IgG_bt2_hg38_rpgc_bin500.bw"
bw_comp="06_ratio/H3K9me3_vs_IgG_log2ratio.bw"

peak_k9="04_peak/H3K9me3-1/H3K9me3-1_peaks.broadPeak"
peak_igg="04_peak/IgG/IgG_peaks.broadPeak"
peak_comp="04_peak/H3K9me3-1_specific.bed"

mkdir 05_peak_against_igg
macs3 callpeak -q 0.1 -t ${sortbam_k9} \
-g hs -n "H3K9me3-1" --broad --broad-cutoff 0.1 \
-c ${sortbam_igg} \
--outdir "05_peak_against_igg/"

mkdir 06_ratio
bigwigCompare -b1 ${bw_k9} -b2 ${bw_igg} --operation log2 -o ${bw_comp}

bedtools intersect -a ${peak_k9} -b ${peak_igg} -f 0.5 -v > ${peak_comp}
