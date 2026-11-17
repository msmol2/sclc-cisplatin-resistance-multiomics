# Processing commands

The analysis ran on the course server, one sample at a time, inside `tmux`, with the shared variables kept in a sourced file (`atac_vars.sh`). The original server scripts are not in this repository. The commands below were rewritten from the tool versions and parameters recorded in the lab report, so options not recorded there (output names, log paths) are illustrative. Paths are shortened. `$S` is the sample name (e.g. `XT_H16_CR_Vor_r1`), and `$T` is the number of threads (16).

## 0. Combined reference

```bash
# prefix fly contigs so they can't collide with human names, then concatenate
sed 's/^>/>dm6_/' dm6.fa > dm6_prefixed.fa
cat GRCh38.primary_assembly.genome.fa dm6_prefixed.fa > fluman.fa
samtools faidx fluman.fa
bowtie2-build --threads $T fluman.fa FLUMAN_GRCh38_primary_dm6_prefixed
```

A TSS BED file (78,696 gene TSSs) and a blacklist restricted to standard chromosomes (910 regions) were made from the GENCODE v50 GTF and ENCODE ENCFF356LFX.

## 1. Trimming

```bash
fastp -i ${S}_R1.fastq.gz -I ${S}_R2.fastq.gz \
      -o trimmed/${S}_R1.fq.gz -O trimmed/${S}_R2.fq.gz \
      --detect_adapter_for_pe --length_required 20 -w $T \
      -j logs/${S}.fastp.json -h logs/${S}.fastp.html
```

## 2. Alignment

```bash
bowtie2 -p $T -x FLUMAN_GRCh38_primary_dm6_prefixed \
        -1 trimmed/${S}_R1.fq.gz -2 trimmed/${S}_R2.fq.gz \
        --no-mixed --no-discordant --maxins 2000 2> logs/${S}.bowtie2.log \
  | samtools sort -@ $T -o bam/${S}.sorted.bam
samtools index bam/${S}.sorted.bam
```

## 3. Spike-in fraction

```bash
samtools idxstats bam/${S}.sorted.bam \
  | awk -v s=$S '$1 ~ /^dm6_/ {fly+=$3} $1 !~ /^dm6_/ && $1 != "*" {hs+=$3}
                 END {printf "%s\t%d\t%d\t%.4f\n", s, hs, fly, fly/(hs+fly)}'
```

## 4. Duplicates and filtering

`fixmate` needs a name-sorted BAM.

```bash
samtools sort -n -@ $T -o tmp.name.bam bam/${S}.sorted.bam
samtools fixmate -m tmp.name.bam tmp.fixmate.bam
samtools sort -@ $T -o tmp.coord.bam tmp.fixmate.bam
samtools markdup -@ $T tmp.coord.bam bam/${S}.markdup.bam
samtools index bam/${S}.markdup.bam

# proper pairs, MAPQ >= 30, no unmapped/secondary/QC-fail/duplicate/supplementary,
# human nuclear chromosomes only
CHROMS=$(samtools idxstats bam/${S}.markdup.bam | cut -f1 | grep -v -E '^dm6_|^chrM$|^\*$')
samtools view -b -f 2 -F 1804 -q 30 bam/${S}.markdup.bam $CHROMS > bam/${S}.filtered.bam
samtools index bam/${S}.filtered.bam
```

## 5. QC

```bash
# fragment sizes
bamPEFragmentSize -b bam/*.filtered.bam --maxFragmentLength 800 \
                  --outRawFragmentLengths qc/fragment_sizes.tab -p $T

# bigWigs and TSS enrichment
bamCoverage -b bam/${S}.filtered.bam -o bw/${S}.bw --binSize 10 --normalizeUsing RPKM -p $T
computeMatrix reference-point -S bw/*.bw -R tss.bed -a 2000 -b 2000 --binSize 10 \
              -o qc/tss_matrix.gz -p $T
plotProfile -m qc/tss_matrix.gz -o qc/TSS_enrichment.png --perGroup

# library complexity
preseq lc_extrap -B -P -o complexity/${S}.lc_extrap.txt bam/${S}.markdup.bam

multiqc -o multiqc .
```

## 6. Peaks and FRiP

```bash
macs3 callpeak -t bam/${S}.filtered.bam -f BAMPE -g hs -n $S --outdir peaks \
               --nomodel --nolambda --keep-dup all

TOTAL=$(samtools view -c bam/${S}.filtered.bam)
INPEAK=$(bedtools intersect -u -a bam/${S}.filtered.bam -b peaks/${S}_peaks.narrowPeak | samtools view -c -)
echo -e "$S\t$TOTAL\t$INPEAK\t$(echo "scale=3; $INPEAK/$TOTAL" | bc)" >> qc/FRiP_scores.txt
```

## 7. Consensus peaks and counts

```bash
cat peaks/*_peaks.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > peaks/consensus_peaks.bed
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
     {print "peak_"NR,$1,$2+1,$3,"+"}' peaks/consensus_peaks.bed > peaks/consensus_peaks.saf

featureCounts -F SAF -a peaks/consensus_peaks.saf -o qc/peak_counts.txt \
              -p --countReadPairs -T $T bam/*.filtered.bam
```

## 8. Differential accessibility (R / edgeR)

```r
library(edgeR)
counts <- read.delim("qc/peak_counts.txt", comment.char = "#", row.names = 1)
counts <- counts[, -(1:5)]
info   <- read.delim("PAGI_ATAC_samplesheet.txt")
colnames(counts) <- info$Sample

condition <- factor(info$Condition)   # default alphabetical levels: CR is the reference
treatment <- factor(info$Treatment, levels = c("Co", "RGFP", "Romi", "Vor"))

y <- DGEList(counts)
y <- y[rowSums(cpm(y) > 1) >= 3, , keep.lib.sizes = FALSE]   # 38,748 peaks
y <- normLibSizes(y)                                           # TMM

design <- model.matrix(~ treatment + condition)
y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
res <- glmQLFTest(fit, coef = "conditionP")   # log2(P/CR); see the note below
tab <- topTags(res, n = Inf)$table
write.csv(tab, "qc/DA_CR_vs_P.csv")

plotMDS(y, col = ifelse(condition == "CR", "steelblue", "tomato"),
        pch = c(15, 16, 17, 18)[treatment])
```

With default factor levels `CR` is the reference, so the tested coefficient is P vs CR: positive log2FC means more open in **parental**. The delivered table was labelled "CR vs P", and the group's DESeq2 cross-check (|r| = 0.95, opposite sign) caught the inversion. The default factor levels are the most likely cause. To get log2(CR/P) directly, use `factor(info$Condition, levels = c("P", "CR"))` and `coef = "conditionCR"`.

`filterByExpr` with default settings kept only 69 peaks at this depth, so the simpler CPM filter above was used instead.
