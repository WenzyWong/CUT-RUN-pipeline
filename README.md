# CUT&RUN-pipeline

This is a pipeline for CUT&RUN analysis.

## For pre-experiment

To be noted: The pre experiment data used in this pipeline didn't contain spike-in and replicates.

Run `cutrun_preexp.py` using snakemake by

```bash
snakemake -s cutrun_preexp.py -p -j ${JOB_NUM} --rerun-incomplete
```

where `${JOB_NUM}` is the number of jobs.

There are some tips that should keep in mind later:

> [!TIP]
> * The sequences in `rule trim` are based on the experimental kit. Therefore, we shall change it according to the kit manual.
> * In `rule bt2_mapping`, using `--very-sensitive`, `--no-mixed`, and `--no-discordant` parameters is due to protocol from [CUT&Tag Data Processing and Analysis Tutorial](https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=9). These parameters would remove background noise as much as possible, thus might generating relatively lower peaks than default `bowtie2` parameters.
> * In `rule bam_to_bw_rpgc`, choosing RPGC as the normalization method used in `bamCoverage` is because:
  > * RPGC normalizes data to a 1x genome coverage, which is important for comparing signal intensities between different samples, especially when sequencing depths vary.
  > * CUT&RUN generates relatively low background noise. It is primarily enriched at true binding sites.
> * In `rule bam_to_bw_rpgc`, setting `--binSize` to 500 is due to relatively large domians of histone modifications. Suitable `binSize` can range from 200 to 1000, according to the researchers' interests.
> * Choosing `macs3` instead of `seacr` to call peaks is due to the reliability of the former method. The `--broad` and `--broad-cutoff 0.1` parameters are requisite for histone modification peak calling.

After the pipeline is finished, run `single_command.sh` by

```bash
bash single_command.sh
```

The commands in this script contain peak-calling treating IgG as input and treating IgG as a control experiment. The former is referenced from [CUT&Tag Data Processing and Analysis Tutorial](https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=31).Personally speaking, I prefer the latter because the experimental processes of generating input in ChIP-seq and generating IgG in CUT&RUN differ a lot. Treating IgG as input would lose many signals with relatively lower peak intensities and cause severe false negative results.

As the former could also provide the most significant differential peaks of treatment sample, I still kept the command to generate a reference peak calling result.

After running these two script, we would get all results like the following structure,

```plain
├── 01_clean
│   ├── H3K9me3-1
│   │   ├── H3K9me3-1_R1.fq.gz
│   │   └── H3K9me3-1_R2.fq.gz
│   └── IgG
│       ├── IgG_R1.fq.gz
│       └── IgG_R2.fq.gz
├── 02_bam
│   ├── H3K9me3-1
│   │   ├── H3K9me3-1_bt2_hg38.sam
│   │   ├── H3K9me3-1_bt2_hg38_sort.bam
│   │   └── H3K9me3-1_bt2_hg38_sort.bam.bai
│   └── IgG
│       ├── IgG_bt2_hg38.sam
│       ├── IgG_bt2_hg38_sort.bam
│       └── IgG_bt2_hg38_sort.bam.bai
├── 03_bw
│   ├── H3K9me3-1
│   │   └── H3K9me3-1_bt2_hg38_rpgc_bin500.bw
│   └── IgG
│       └── IgG_bt2_hg38_rpgc_bin500.bw
├── 04_peak
│   ├── H3K9me3-1
│   │   ├── H3K9me3-1_model.r
│   │   ├── H3K9me3-1_peaks.broadPeak
│   │   ├── H3K9me3-1_peaks.gappedPeak
│   │   └── H3K9me3-1_peaks.xls
│   ├── H3K9me3-1_specific.bed
│   └── IgG
│       ├── IgG_model.r
│       ├── IgG_peaks.broadPeak
│       ├── IgG_peaks.gappedPeak
│       └── IgG_peaks.xls
├── 05_peak_against_igg
│   ├── H3K9me3-1_model.r
│   ├── H3K9me3-1_peaks.broadPeak
│   ├── H3K9me3-1_peaks.gappedPeak
│   └── H3K9me3-1_peaks.xls
├── 06_ratio
│   ├── H3K9me3_ratio_peaks.bed
│   ├── H3K9me3_vs_IgG_log2ratio.bw
└──logs
    ├── 01_clean
    │   ├── H3K9me3-1_cutadapt.log
    │   └── IgG_cutadapt.log
    ├── 02_bam
    │   ├── H3K9me3-1_bt2_hg38.log
    │   ├── H3K9me3-1_bt2_hg38_sort.log
    │   ├── IgG_bt2_hg38.log
    │   └── IgG_bt2_hg38_sort.log
    ├── 03_bw
    │   ├── H3K9me3-1_bw_rpgc_bin500.log
    │   └── IgG_bw_rpgc_bin500.log
    ├── 04_peak
    │   ├── H3K9me3-1_peaks.log
    │   └── IgG_peaks.log
    └── 05_peak_against_igg
        └── H3K9me3-1_peaks.log
```

## For experiment with spike-in

To be noted: Spike-in and replicates are included in this pipeline.

Run `cutrun_spikein.py` using snakemake by

```bash
snakemake -s cutrun_spikein.py -p -j ${JOB_NUM} --rerun-incomplete
```

where `${JOB_NUM}` is the number of jobs.

Tips can be found in the [previous section](#for-pre-experiment).
