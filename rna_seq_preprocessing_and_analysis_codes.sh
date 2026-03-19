fastqc -o fastq_raw -t 8 *.fastq.gz
# -o - direktorija, į kuria saugos html failus (fastq_raw) šita direktorija turi būti sukurta iš anksto
mkdir fastq_trimmed
# fastq_trimmed yra direktorijoje raw_data ir trim_galore naudojami originlūs fastq.gz failai
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647688_1.fastq.gz SRR11647688_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647689_1.fastq.gz SRR11647689_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647690_1.fastq.gz SRR11647690_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647699_1.fastq.gz SRR11647699_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647700_1.fastq.gz SRR11647700_2.fastq.gz
trim_galore -o fastq_trimmed -j 8 --quality 20 --length 20 --paired SRR11647701_1.fastq.gz SRR11647701_2.fastq.gz
# -j 8 (galima naudoti ir 6, kuo didesnis skaičius, tuo greičiau) pagreitina veiksmą
cd fastq_trimmed
mkdir fastqc_trimmed
fastqc -o fastqc_trimmed -t 8 *.fq.gz
cd fastqc_trimmed
multiqc .
# multiqc trimmintiems html failams
cd fastq_raw
multiqc .
# multiqc pradiniams html 
# references direktorijoje:
gunzip GRCh38.primary_assembly.genome.fa.gz
# Homo sapiens ir dideliems genomams padeda išvengti klaidų (tada nebe .fa.gz, bet lieka .fa)
hisat2-build GRCh38.primary_assembly.genome.fa indexed_ref
# reikėtų pridėti -p 6 arba 8
cd ~/HW1/raw_data
mkdir mapped
cd mapped
hisat2 -p 6 -q -x ~/HW1/references/indexed_ref -1 ~/HW1/raw_data/fastq_trimmed/SRR11647688_1_val_1.fq.gz -2 ~/HW1/raw_data/fastq_trimmed/SRR11647688_2_val_2.fq.gz -S  SRR11647688.sam
hisat2 -p 6 -q -x ~/HW1/references/indexed_ref -1 ~/HW1/raw_data/fastq_trimmed/SRR11647689_1_val_1.fq.gz -2 ~/HW1/raw_data/fastq_trimmed/SRR11647689_2_val_2.fq.gz -S  SRR11647689.sam
hisat2 -p 6 -q -x ~/HW1/references/indexed_ref -1 ~/HW1/raw_data/fastq_trimmed/SRR11647690_1_val_1.fq.gz -2 ~/HW1/raw_data/fastq_trimmed/SRR11647690_2_val_2.fq.gz -S  SRR11647690.sam
hisat2 -p 6 -q -x ~/HW1/references/indexed_ref -1 ~/HW1/raw_data/fastq_trimmed/SRR11647699_1_val_1.fq.gz -2 ~/HW1/raw_data/fastq_trimmed/SRR11647699_2_val_2.fq.gz -S  SRR11647699.sam
hisat2 -p 6 -q -x ~/HW1/references/indexed_ref -1 ~/HW1/raw_data/fastq_trimmed/SRR11647700_1_val_1.fq.gz -2 ~/HW1/raw_data/fastq_trimmed/SRR11647700_2_val_2.fq.gz -S  SRR11647700.sam
hisat2 -p 6 -q -x ~/HW1/references/indexed_ref -1 ~/HW1/raw_data/fastq_trimmed/SRR11647701_1_val_1.fq.gz -2 ~/HW1/raw_data/fastq_trimmed/SRR11647701_2_val_2.fq.gz -S  SRR11647701.sam
# galima pridėti --no-unal --no-discordant --no-mixed (nebūtina)
# Number of reads used for mapping:
samtools flagstat SRR11647688.sam
samtools flagstat SRR11647689.sam
samtools flagstat SRR11647690.sam
samtools flagstat SRR11647699.sam
samtools flagstat SRR11647700.sam
samtools flagstat SRR11647701.sam
mkdir bam_files
samtools view -bS SRR11647688.sam > ~/HW1/raw_data/mapped/bam_files/SRR11647688.bam
samtools view -bS SRR11647689.sam > ~/HW1/raw_data/mapped/bam_files/SRR11647689.bam
samtools view -bS SRR11647690.sam > ~/HW1/raw_data/mapped/bam_files/SRR11647690.bam
samtools view -bS SRR11647699.sam > ~/HW1/raw_data/mapped/bam_files/SRR11647699.bam
samtools view -bS SRR11647700.sam > ~/HW1/raw_data/mapped/bam_files/SRR11647700.bam
samtools view -bS SRR11647701.sam > ~/HW1/raw_data/mapped/bam_files/SRR11647701.bam
cd bam_files
mkdir duplicates
mkdir correlation
mkdir PCA
mkdir coverage
mkdir gene_coverage
mkdir inner_distance
mkdir clipping_profile
mkdir annotated_junctions
# duplikatai:
samtools collate -@ 4 -O -u SRR11647688.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/raw_data/mapped/bam_files/duplicates/SRR11647688.bam
samtools collate -@ 4 -O -u SRR11647689.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/raw_data/mapped/bam_files/duplicates/SRR11647689.bam
samtools collate -@ 4 -O -u SRR11647690.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/raw_data/mapped/bam_files/duplicates/SRR11647690.bam
samtools collate -@ 4 -O -u SRR11647699.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/raw_data/mapped/bam_files/duplicates/SRR11647699.bam
samtools collate -@ 4 -O -u SRR11647700.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/raw_data/mapped/bam_files/duplicates/SRR11647700.bam
samtools collate -@ 4 -O -u SRR11647701.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ~/HW1/raw_data/mapped/bam_files/duplicates/SRR11647701.bam
# duplikatų kiekio patikrinimas:
cd duplicates
samtools view -c -f 1024 SRR11647688.bam
samtools view -c -f 1024 SRR11647689.bam
samtools view -c -f 1024 SRR11647690.bam
samtools view -c -f 1024 SRR11647699.bam
samtools view -c -f 1024 SRR11647700.bam
samtools view -c -f 1024 SRR11647701.bam
# Gene body coverage:
cd ~/HW1/raw_data/mapped/bam_files
samtools sort -@ 6 SRR11647688.bam -o ~/HW1/raw_data/mapped/bam_files/gene_coverage/SRR11647688.sorted.bam 
samtools sort -@ 6 SRR11647689.bam -o ~/HW1/raw_data/mapped/bam_files/gene_coverage/SRR11647689.sorted.bam
samtools sort -@ 6 SRR11647690.bam -o ~/HW1/raw_data/mapped/bam_files/gene_coverage/SRR11647690.sorted.bam
samtools sort -@ 6 SRR11647699.bam -o ~/HW1/raw_data/mapped/bam_files/gene_coverage/SRR11647699.sorted.bam
samtools sort -@ 6 SRR11647700.bam -o ~/HW1/raw_data/mapped/bam_files/gene_coverage/SRR11647700.sorted.bam
samtools sort -@ 6 SRR11647701.bam -o ~/HW1/raw_data/mapped/bam_files/gene_coverage/SRR11647701.sorted.bam
# samtools naudoja -@ vietoj -p, kad būtų greičiau
# indeksavimas:
cd gene_coverage
samtools index SRR11647688.sorted.bam
samtools index SRR11647689.sorted.bam
samtools index SRR11647690.sorted.bam
samtools index SRR11647699.sorted.bam
samtools index SRR11647700.sorted.bam
samtools index SRR11647701.sorted.bam
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE_V47.bed.gz/download
mv download hg38_GENCODE_V47.bed.gz
gunzip hg38_GENCODE_V47.bed.gz
# Grafiko padarymas:
conda activate rseqc_env
geneBody_coverage.py -r hg38_GENCODE_V47.bed -i *.sorted.bam -o prefix 
# Inner distance: 
# Inner distance parodo atstumų tarp paired-end reads pasiskirstymą
cd ~/HW1/raw_data/mapped/bam_files/gene_coverage
for bam in ~/HW1/raw_data/mapped/bam_files/gene_coverage/*sorted.bam; do
    samtools view -f 2 "$bam" | awk '{if($9>0) print $9}' >> ~/HW1/raw_data/mapped/bam_files/inner_distance/all_insert_sizes.txt
done
# Per R:
cd ~/HW1/raw_data/mapped/bam_files/inner_distance
R
data <- scan("~/HW1/raw_data/mapped/bam_files/inner_distance/all_insert_sizes.txt")
data <- data[data > 0 & data < 1000]
png("~/HW1/raw_data/mapped/bam_files/inner_distance/inner_distance_distribution.png",
    width=900, height=600)
hist(data,
     breaks=100,
     main="Inner Distance Distribution",
     xlab="Insert Size (bp)",
     ylab="Frequency")
dev.off()
# Clipping profile:
#!/bin/bash
IN=~/HW1/raw_data/mapped/bam_files/gene_coverage
OUT=~/HW1/raw_data/mapped/bam_files/clipping_profile
cd "$OUT"
conda activate rseqc_env
for b in "$IN"/*.sorted.bam; do
    n=$(basename "$b" .sorted.bam)
    clipping_profile.py -i "$b" -o "$OUT/${n}_clipping" -s PE
done
cd ~/HW1/raw_data/mapped/bam_files/clipping_profile
for f in *.pdf; do
    pdftoppm "$f" "${f%.pdf}" -png
done
# Annotated junctions:
cd ~/HW1/raw_data/mapped/bam_files
#!/bin/bash
IN=~/HW1/raw_data/mapped/bam_files/gene_coverage
OUT=~/HW1/raw_data/mapped/bam_files/annotated_junctions
REF=$IN/hg38_GENCODE_V47.bed
for b in "$IN"/*.sorted.bam; do
    n=$(basename "$b" .sorted.bam)
    junction_annotation.py -i "$b" -r "$REF" -o "$OUT/${n}_junction"
done
cat "$OUT"/*_junction.bed | sort -k1,1 -k2,2n > "$OUT/all.bed"
b=$(ls "$IN"/*.sorted.bam | head -n1)
junction_annotation.py -i "$b" -r "$REF" -o "$OUT/combined"
cp "$OUT/all.bed" "$OUT/combined.junction.bed"
junction_annotation.py -i "$b" -r "$REF" -o "$OUT/combined_final"
for f in *.pdf; do
    pdftoppm "$f" "${f%.pdf}" -png
done
# Koreliacija:
cd correlation
multiBamSummary bins -p 6 --bamfiles ~/HW1/raw_data/mapped/bam_files/gene_coverage/*.sorted.bam -o correlation_matrix.npz 
plotCorrelation -in correlation_matrix.npz -c spearman -p heatmap -o plot.png 
# Coverage:
cd coverage
samtools coverage ~/HW1/raw_data/mapped/bam_files/SRR11647689.bam
samtools coverage ~/HW1/raw_data/mapped/bam_files/SRR11647689.bam
samtools coverage ~/HW1/raw_data/mapped/bam_files/SRR11647690.bam
samtools coverage ~/HW1/raw_data/mapped/bam_files/SRR11647699.bam
samtools coverage ~/HW1/raw_data/mapped/bam_files/SRR11647700.bam
samtools coverage ~/HW1/raw_data/mapped/bam_files/SRR11647701.bam
# PCA:
cd PCA
plotPCA -in ~/HW1/raw_data/mapped/bam_files/correlation/correlation_matrix.npz -o pca.png
# Read counting:
cd ~/HW1/raw_data/mapped/bam_files
mkdir read_counts
cd read_counts
featureCounts -a ~/HW1/references/gencode.v49.primary_assembly.basic.annotation.gtf.gz -o counts.txt -t exon -g gene_id -p --countReadPairs -s 0 ~/HW1/raw_data/mapped/bam_files/SRR*.bam
# Kad atsidaryti HTML failus:
# cd į ten kur yra HTML failai
python -m http.server 8000
# ir tada galima atsidaryti per browser; kad uždaryti Ctrl C
# Chromosomų kiekio patikrinimas reference faile:
samtools faidx GRCh38.primary_assembly.genome.fa
cut -f1,2 GRCh38.primary_assembly.genome.fa.fai
# Determine the library type (stranded or unstranded; if stranded, indicate which strand is first):
infer_experiment.py -r ~/HW1/raw_data/mapped/bam_files/gene_coverage/hg38_GENCODE_V47.bed -i ~/HW1/raw_data/mapped/bam_files/gene_coverage/SRR11647688.sorted.bam
# This is PairEnd Data
# Fraction of reads failed to determine: 0.3124
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.3487
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.3388
# Create a plot (using any tool of your choice) showing the mapping rate for each sample:
# Per R:
R
library(ggplot2)
library(tidyverse)
library(readr)
data_table <- data.frame(
  Sample = c("SRR11647688", "SRR11647689", "SRR11647690",
             "SRR11647699", "SRR11647700", "SRR11647701"),
  Unique = c(90.96, 89.44, 89.52, 91.32, 83.46, 91.55),
  NonUnique = c(4.15, 4.37, 3.34, 3.01, 10.66, 3.17)
)
data_long <- data_table %>%
  pivot_longer(cols = c("Unique", "NonUnique"),
               names_to = "MappingType",
               values_to = "Percentage")
data_long$MappingType <- factor(data_long$MappingType, levels = c("NonUnique", "Unique"))
data_long <- data_long %>%
  group_by(Sample) %>%
  arrange(MappingType) %>%
  mutate(pos = cumsum(Percentage) - Percentage / 2)                    
mapping_rate <- ggplot(data_long, aes(x = Sample, y = Percentage, fill = MappingType)) +
  geom_bar(stat = "identity") +
  labs(title = "Mapping Rate",
       x = "Sample",
       y = "Percentage of Reads (%)",
       fill = "Mapping type") +
  scale_fill_manual(
    values = c("Unique" = "lightblue", "NonUnique" = "#C080FF"),
    labels = c("Non-unique", "Unique")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white", color = NA)
    )
ggsave("Mapping_rate.png", plot = mapping_rate, width = 8, height = 6, dpi = 300, bg = "white")
# Create a stacked bar plot showing feature assignment rates (from featureCounts):
# Per R:
R
library(ggplot2)
library(tidyverse)
library(readr)
counts_txt <- read_table("counts.txt.summary")
colnames(counts_txt) <- c("Status", "SRR11647688", "SRR11647689", "SRR11647690", "SRR11647699", "SRR11647700", "SRR11647701")
counts_table_filtered <- counts_txt[rowSums(counts_txt[,-1]) > 0, ]
counts_txt_long <- pivot_longer(counts_table_filtered, -Status,
                        names_to="Sample",
                        values_to="Count")
feature_rates <- ggplot(counts_txt_long, aes(x = Sample, y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Read Assignment Summary",
       y = "Read Counts",
       x = "Sample")
ggsave("Read_assignment_summary.png", plot = feature_rates, width = 8, height = 6, dpi = 300, bg = "White")
# Analizė per R:
# Data statistical analysis
# Prieš analizę per R:
cd ~/HW1/raw_data/mapped/bam_files/read_counts
zcat ~/HW1/references/gencode.v49.primary_assembly.basic.annotation.gtf.gz \ | awk 'BEGIN{FS="\t"; OFS="\t"} $3=="gene" {
    gene_id=""; gene_type="";
    if (match($9, /gene_id "[^"]+"/)) {
        gene_id=substr($9, RSTART+9, RLENGTH-10)
    }
    if (match($9, /gene_type "[^"]+"/)) {
        gene_type=substr($9, RSTART+11, RLENGTH-12)
    }
    print gene_id, gene_type
}' > gene_info.tsv
R
# Lentelių parsisiuntimui:
library(readr)
setwd("~/HW1/raw_data/mapped/bam_files/read_counts")
counts_raw <- read_table("counts.txt", 
    skip = 1)
gene_types <- read_table("gene_info.tsv", 
    col_names = FALSE)
colnames(gene_types) <- c("Gene_ID", "Gene_type")
counts_clean <- counts_raw[, -c(2:6)]
colnames(counts_clean) <- c("Gene_ID", "SRR11647688", "SRR11647689", "SRR11647690", "SRR11647699", "SRR11647700", "SRR11647701")
# Column_data reikia susikurti lentelę:
column_data <- data.frame(
  Sample = c("SRR11647688", "SRR11647689", "SRR11647690",
             "SRR11647699", "SRR11647700", "SRR11647701"),
  Type = c("Non-tumor", "Non-tumor", "Non-tumor",
                "Tumor", "Tumor", "Tumor")
)
print(column_data) 
write.table(column_data,
            file = "column_data.txt",
            sep = " ",
            row.names = FALSE,
            quote = FALSE)
print(column_data)
library(readxl)
coldata <- read.table("column_data.txt", header=TRUE)
rownames(coldata) <- coldata$Sample
coldata$Type <- as.factor(coldata$Type)
print(coldata)
# Duomenų sorting ir matricos padarymas:
combined_table <- merge(gene_types, counts_clean, by = "Gene_ID")
protein_coding <- combined_table[combined_table$Gene_type == "protein_coding", ]
counts_protein_coding <- protein_coding[, -2]
rownames(counts_protein_coding) <- counts_protein_coding$Gene_ID
counts_protein_coding <- counts_protein_coding[, -1]
# Differential Expression su DESeq2 paketu:
# install.packages("BiocManager", lib="~/R/library")
library(BiocManager)
 .libPaths("~/R/library")
BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = counts_protein_coding,
  colData = coldata,
  design = ~ Type
)
dds$Type <- relevel(dds$Type, ref = "Non-tumor")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes remaining after low-count filter:", nrow(dds), "\n")
# Genes remaining after low-count filter: 16554
# Pasiliekamos tik eilutės, kuriose suma eiluteje yra daugiau nei 10 (todėl, nes nereikšmingas padidejimas klaidina)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Type","Tumor","Non-tumor"), alpha = 0.05)
summary(res)
# out of 16227 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 791, 4.8%
# LFC < 0 (down)     : 1179, 7.1%
# outliers [1]       : 174, 1.1%
# low counts [2]     : 321, 1.9%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
.libPaths("~/R/library")
BiocManager::install("apeglm", update = FALSE, ask = FALSE)
library(apeglm)
res_shrunk <- lfcShrink(dds, coef = "Type_Tumor_vs_Non.tumor", type = "apeglm")
# using 'apeglm' for LFC shrinkage. If used in published research, please cite:
   # Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
   # sequence count data: removing the noise and preserving large differences.
   # Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
# Lentelės reikalingos fold-change:
library(dplyr)
library(tibble)
res_df <- res_shrunk |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  arrange(padj) |>
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"  # Not Significant
    ),
    significance = factor(significance, levels = c("Up", "Down", "NS"))
  )
table(res_df$significance)
#  Up  Down    NS 
# 713  1116 14725 
# Atrenkami pirmi 10 turintys didžiausią fold-change:
cat("\nTop 10 upregulated genes:\n")
res_df %>% filter(significance == "Up") %>% head(10) %>%
  dplyr::select(gene, log2FoldChange, padj) %>% print()
#   gene log2FoldChange         padj
# 1  ENSG00000060718.23      10.072689 8.056855e-46
# 2  ENSG00000099953.11       7.937720 2.620550e-35
# 3  ENSG00000123500.10       8.576835 4.324182e-33
# 4  ENSG00000151388.11       5.736830 1.843922e-31
# 5  ENSG00000137745.14       7.300126 4.548444e-28
# 6  ENSG00000133110.16       6.171963 4.726416e-25
# 7  ENSG00000106483.13       5.943427 2.858737e-24
# 8  ENSG00000164694.18       5.564844 7.851905e-22
# 9  ENSG00000108821.15       5.861209 3.067107e-21
# 10 ENSG00000148848.15       5.553848 3.842370e-20
# Be fold change:
deg_no_fc <- res_df %>%
  filter(padj < 0.05)
write.csv(deg_no_fc, "DEG_no_fold_change.csv", row.names = FALSE)
# Fold change:
deg_fc <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg_fc, "DEG_with_fold_change.csv", row.names = FALSE)
# Grafikai:
# 1. Identify total number of reads per sample and number of detected genes:
reads_per_sample <- colSums(counts_clean[, -1])
genes_detected <- colSums(counts_clean[, -1] > 0)
qc_table <- data.frame(
  Sample = colnames(counts_clean[, -1]),
  Total_reads = reads_per_sample,
  Detected_genes = genes_detected
)
library(tidyr)
qc_table_long <- qc_table |>
  pivot_longer(cols = c(Total_reads, Detected_genes),
               names_to = "Metric",
               values_to = "Value") |>
  mutate(log_value = log2(Value + 1))
library(ggplot2)
# 1 grafikas:
plot1 <- qc_table_long %>%
  filter(Metric == "Total_reads") %>%
  ggplot(aes(x = Sample, y = Value)) +
  geom_bar(stat = "identity", fill = "#BA00C4") +
  theme_minimal() +
  labs(
    title = "Total Reads per Sample",
    x = "Sample",
    y = "Total Reads"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )
ggsave("total_reads_plot.png", plot = plot1, width = 8, height = 6, dpi = 300, bg = "white")
# 2 grafikas:
plot2 <- qc_table_long %>%
  filter(Metric == "Detected_genes") %>%
  ggplot(aes(x = Sample, y = Value)) +
  geom_bar(stat = "identity", fill = "#6D99F8") +
  theme_minimal() +
  labs(
    title = "Detected Genes per Sample",
    x = "Sample",
    y = "Number of Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )
ggsave("detected_genes_plot.png", plot = plot2, width = 8, height = 6, dpi = 300, bg = "white")
# 2.-3. Make a correlation plot that shows pearson correlation (using gene level data) between individual samples and provide a PCA plot:
dds_all <- DESeqDataSetFromMatrix(
  countData = counts_protein_coding[rowSums(counts_protein_coding) >= 10, ],
  colData = coldata,
  design = ~ Type
)
vst_all <- vst(dds_all, blind = FALSE)
norm_counts <- assay(vst_all)
write.csv(norm_counts, "normalized_counts_VST.csv")
# 3 grafikas PCA plot:
library(ggrepel)
pcaData <- plotPCA(vst_all, intgroup = "Type", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot3 <- ggplot(pcaData, aes(PC1, PC2, color = Type, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(color = "Type") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
    )
ggsave("PCA.png", plot = plot3, width = 8, height = 6, dpi = 300, bg = "white")
# 4 grafikas heatmap:
.libPaths("~/R/library")
install.packages("pheatmap", repos = "https://cloud.r-project.org", lib="~/R/library")
library(pheatmap)
cor_mat  <- cor(norm_counts, method = "pearson")
# Annotation for the heatmap
anno_df <- coldata |>
  dplyr::select(Type)
anno_colors <- list(
  Type = c(
    Tumor                 = "#F74062",
    'Non-tumor'                 = "#40F76B"
  )
)
png("correlation.png", width = 800, height = 600, res = 150)
pheatmap(
  cor_mat,
  annotation_col = anno_df,
  annotation_row = anno_df,
  annotation_colors = anno_colors,
  color = colorRampPalette(c("white", "#55036B"))(50),
  border_color = NA,
  main = "Sample-to-Sample Pearson Correlation",
  fontsize = 10
)
dev.off()
# 4. Compare raw and normalized counts visually:
# 5 grafikas boxplot:
png("boxplot_raw_counts.png", width = 800, height = 600)
par(mar = c(7, 4, 3, 2))
plot4 <- boxplot(log2(counts_protein_coding + 1),
        las = 2,
        col = "#C8A2C8",
        main = "Raw counts",
        ylab = "log2(count + 1)")
dev.off()
# 6 grafikas boxplot normalized:
png("vst_normalized_counts.png", widt = 800, height = 600)
par(mar = c(7, 4, 3, 2))
boxplot(assay(vst_all),
        las = 2,
        col = "lightblue",
        main = "VST normalized counts",
        ylab = "Normalized expression")
dev.off()
# 5. Provide a volcano plot:
# 7 grafikas volcano plot:
volcano_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(negLog10Padj = -log10(padj))
top_labels <- volcano_df %>%
  filter(significance != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 10)
plot4 <- ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10Padj, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_text_repel(
    data = top_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = 100,
    box.padding = 0.3
  ) +
  scale_color_manual(values = c("Up" = "#AA00C4", "Down" = "#0089C4", "NS" = "grey70")) +
  labs(
    title = "Volcano Plot — tumor vs non-tumor",
    subtitle = "padj < 0.05 and |log2FC| > 1 threshold shown as dashed lines",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "darkgray", fill = NA, linewidth = 1),
    plot.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )
ggsave("Volcano_plot.png", plot = plot4, width = 8, height = 6, dpi = 300, bg = "white")
# 6. Create/provide an MA plot:
# 8 grafikas MA plot:
png("MA_plot.png", width = 800, height = 600)
DESeq2::plotMA(res_shrunk, ylim = c(-6, 6), main = "MA plot — tumor vs non-tumor (apeglm shrinkage)")
dev.off()
# 9 grafikas MA plot shrunk:
png("MA_plot_not_shrunk.png", width = 800, height = 600)
DESeq2::plotMA(res, ylim = c(-6, 6), main = "MA plot — tumor vs non-tumor (apeglm shrinkage)")
dev.off()
# 7. Create/provide heatmap:
# 10 grafikas heatmap:
# Take top 50 DE genes (by padj) for the Listerine vs Untreated comparison
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)
# Subset the VST matrix to Listerine + untreated samples only, top genes
samples_tumor_nontumor <- coldata %>%
  filter(Type %in% c("Tumor", "Non-tumor")) %>%
  pull(Sample)
vst_subset <- norm_counts[top_genes, samples_tumor_nontumor]
# Z-score normalize by row (gene) so we see relative expression pattern
vst_z <- t(scale(t(vst_subset)))
# Annotation
anno_subset <- coldata %>%
  filter(Sample %in% samples_tumor_nontumor) %>%
  dplyr::select(Type)
library(RColorBrewer)
png("heatmap.png", width = 800, height = 600, res = 150)
pheatmap(
  vst_z,
  annotation_col  = anno_subset,
  annotation_colors = list(treatment = anno_colors$Type[c("Tumor", "Non-tumor")]),
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  border_color = NA,
  main = "Top 50 DE Genes - Tumor vs Non-tumor (Z-scored VST)"
)
dev.off()
# 8. Compare raw p-value versus adjusted p-values:
# 11 grafikas:
png("p.values.png", width = 800, height = 600, res = 150)
par(mfrow = c(1, 2))
hist(res_df$pvalue,
     breaks = 50,
     col = "#C8A2C8",
     main = "Raw p-values",
     xlab = "p-value")
hist(res_df$padj,
     breaks = 50,
     col = "lightblue",
     main = "Adjusted p-values",
     xlab = "adjusted p-value")
dev.off()
sum(res_df$pvalue < 0.05, na.rm = TRUE)
# [1] 3636
sum(res_df$padj < 0.05, na.rm = TRUE)
# [1] 1986
# Data biological analysis:
# 1. Pick TOP 3 upregulated genes and TOP 3 downregulated genes:
top3_upregulated <- res_df %>%
  filter(significance == "Up") %>%
  arrange(padj) %>%
  slice_head(n = 3)
top3_downregulated <- res_df %>%
  filter(significance == "Down") %>%
  arrange(padj) %>%
  slice_head(n = 3)
top3_upregulated <- merge(top3_upregulated, protein_coding, by = 1)
top3_downregulated <- merge(top3_downregulated, protein_coding, by = 1)
print(top3_upregulated)
#  gene baseMean log2FoldChange     lfcSE       pvalue
# 1 ENSG00000060718.23 2953.583      10.072689 0.6897191 2.559682e-49
# 2 ENSG00000099953.11 1638.109       7.937720 0.6213099 1.665110e-38
# 3 ENSG00000123500.10 1967.347       8.576835 0.6942647 3.022367e-36
#           padj significance      Gene_type SRR11647688 SRR11647689 SRR11647690
# 1 8.056855e-46           Up protein_coding           5           8           2
# 2 2.620550e-35           Up protein_coding           2          10          25
# 3 4.324182e-33           Up protein_coding           7          20           1
#   SRR11647699 SRR11647700 SRR11647701
# 1       10859        4569        2645
# 2        3334        4033        2532
# 3        6375        3046        2638
print(top3_downregulated)
# gene  baseMean log2FoldChange     lfcSE       pvalue
# 1 ENSG00000096006.12  9718.325      -12.20668 0.7587199 6.693132e-59
# 2 ENSG00000117983.17 27424.642      -13.82778 0.8618146 1.590169e-58
# 3  ENSG00000204544.5 13579.332      -11.21142 0.7254577 6.926010e-55
#           padj significance      Gene_type SRR11647688 SRR11647689 SRR11647690
# 1 1.053365e-54         Down protein_coding       10226       37230       10100
# 2 1.251304e-54         Down protein_coding       16965      118811       27966
# 3 3.633385e-51         Down protein_coding       33720       32486       12111
#   SRR11647699 SRR11647700 SRR11647701
# 1           1           2           9
# 2           3           0           8
# 3           5           2          27
# 2. Perform GO over-representation analysis:
sig_de <- res_df |>
  filter(padj < 0.05 & abs(log2FoldChange) > 1)
library(stringr)
sig_de$gene <- str_remove(sig_de$gene, "\\..*")
sig_de_genes <- sig_de$gene
sig_de_genes <- unique(sig_de_genes[!is.na(sig_de_genes)])
library(clusterProfiler)
library(org.Hs.eg.db)
ego_all <- enrichGO(
  gene          = sig_de_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
# 12-14 grafikai pakitusiems keliams:
# 12 dot plot:
plot_GO <- dotplot(ego_all, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free")
ggsave("GO_dotplot_faceted.png", plot = plot_GO, width = 8, height = 10, dpi = 300)
# 13 bar plot:
plot_GO_bar <- barplot(ego_all, showCategory = 20)
ggsave("GO_barplot.png", plot = plot_GO_bar, width = 8, height = 10, dpi = 300)
# 14 emap:
library(enrichplot)
plot_emap <- emapplot(pairwise_termsim(ego_all), showCategory = 20)
ggsave("GO_emapplot.png", plot = plot_emap, width = 10, height = 8, dpi = 300)
# 3. Using your list of genes, perform GSEA analysis over GO and MSIGdb:
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
.libPaths("~/R/library")
install.packages("msigdbr", repos = "https://cloud.r-project.org", lib="~/R/library")
library(msigdbr)
library(BiocParallel)
library(enrichplot)
BiocParallel::register(BiocParallel::SerialParam())
# geneList paruošimas GSEA analizei:
gsea_geneList <- res_df %>%
  dplyr::mutate(gene = sub("\\..*", "", gene)) %>%   # pašalina ENSEMBL versijas po taško
  dplyr::filter(!is.na(gene), !is.na(log2FoldChange)) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FoldChange) %>%
  tibble::deframe() %>%
  sort(decreasing = TRUE)
# papildomi patikrinimai:
gsea_geneList <- gsea_geneList[!is.na(names(gsea_geneList))]
gsea_geneList <- gsea_geneList[!duplicated(names(gsea_geneList))]
gsea_geneList <- sort(gsea_geneList, decreasing = TRUE)
# GO GSEA:
gsea_go <- clusterProfiler::gseGO(
  geneList = gsea_geneList,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  verbose = FALSE
)
# MSigDB Hallmark:
msig_h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  dplyr::filter(!is.na(ensembl_gene)) %>%
  dplyr::distinct()
gsea_hallmark <- clusterProfiler::GSEA(
  geneList = gsea_geneList,
  TERM2GENE = msig_h,
  pvalueCutoff = 0.05,
  verbose = FALSE
)
# 1 gsea grafikas:
plot_gsea <- gseaplot2(gsea_go, geneSetID = 1:4)
ggsave("GSEA_enrichment_plot.png", plot = plot_gsea, width = 8, height = 6, dpi = 300)
# 2 gsea grafikas:
plot_gsea2 <- gseaplot2(gsea_hallmark, geneSetID = 1:4)
ggsave("GSEA_hallmark_plot.png", plot = plot_gsea2, width = 8, height = 6, dpi = 300)
# 3 ridgeplot:
.libPaths("~/R/library")
install.packages("ggridges", repos = "https://cloud.r-project.org", lib="~/R/library")
library(ggridges)
plot_ridge <- ridgeplot(gsea_go, showCategory = 10)
ggsave("Ridgeplot_gsea_go.png", plot = plot_ridge, width = 8, height = 6, dpi = 300)
# 4 ridgeplot
plot_ridge2 <- ridgeplot(gsea_hallmark, showCategory = 10)
ggsave("Ridgeplot_hallmark.png", plot = plot_ridge2, width = 8, height = 6, dpi = 300)
# 5 dot plot:
plot_dot <- dotplot(gsea_go, showCategory = 10)
ggsave("Dotplot_gsea_go.png", plot = plot_dot, width = 8, height = 6, dpi = 300)
# 6 dot plot:
plot_dot2 <- dotplot(gsea_hallmark, showCategory = 10)
ggsave("Dotplot_gsea_hallmark.png", plot = plot_dot2, width = 8, height = 6, dpi = 300)
# Hallmark:
gsea_hallmark <- pairwise_termsim(gsea_hallmark)
# enrichment network plot:
p_hallmark <- emapplot(
  gsea_hallmark,
  showCategory = 20
)
ggsave("Hallmark_enrichemnt_network_plot.png", plot = p_hallmark, width = 8, height = 6, dpi = 300)
# 4. Using aPEAR package, make a map visualisation of your GSEA results:
.libPaths("~/R/library")
library(aPEAR)
# GO GSEA map:
p_go <- enrichmentNetwork(
  gsea_go@result,
  drawEllipses = TRUE,
  fontSize = 4,
  repel = TRUE
)
ggsave("Enrichment_network.png", plot = p_go, width = 8, height = 6, dpi = 300)
p_hallmark <- enrichmentNetwork(
  gsea_hallmark@result,
  drawEllipses = TRUE,
  fontSize = 3,
)
ggsave("Hallmark_enrichment_network.png", plot = p_hallmark, width = 8, height = 6, dpi = 300)
# Error in findClusters(sim, method = clustMethod, nameMethod = clustNameMethod,  : 
  # No clusters found.
# Hallmark GSEA results don't have enough enriched terms to form a network and the plot cannot be made.