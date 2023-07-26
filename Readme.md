# requirements:

## Fastq files (also Fasta)

Fastq (and Fasta) files are text files that contain the individual reads of a sequencing run, in combination with their quality scores. They look like this:

```
@A00643:174:HYTFTDRXX:1:2101:1289:1000 1:N:0:GCCAAT
NAAGCTGTAGTTCCACTTTTCAAACACCAAAAAGGGGGGAAAGGGGGAGTTGTCTTCGCTGTCAAATAGTTACATTGTTAATTCCATAATAGCATTTACAA
+
#FFFFFFFFFFFF:FFFFFFFFF,FFFFFFFF,FFFFFFF,FFFFFFFFFFFF,,:,F:,FFFFF,FFFFFFFFFFFFF:FFF::FFF,F:,F:FFFF,F,
```

1.  sequence identifier (name); always starts with `@` or `>`
2.  sequence (A,C,T,G and N)
3.  seqerator between lines
4.  base quality score as letters

They are often compressed, ending in ".gz". These can be extracted using `gunzip`.
In the case of paired end reads, the file names often contain "R1" and "R2" indicating the mate files.

## Reference genome

The reads have to be mapped against something. The reference genome can be downloaded for example from *ENSEMBL*. Here is a link for the human *GRCh38* genome:
[ENSEMBL: homo sapiens](http://www.ensembl.org/Homo_sapiens/Info/Index)
It contains different version.
First Option:

- dna_sm - Repeats soft-masked (converts repeat nucleotides to lowercase)
- dna_rm - Repeats masked (converts repeats to to N's)
- dna - No masking

second option:

- .toplevel - Includes haplotype information; common variants of the reference
- .primary_assembly - Single reference base per position

I recommend primary_assembly, soft-masked or unmasked. The haplotype information can cause incorrect multi-mapping. But software also has heuristics to discard alignments in repeats, so masking is not necessary.

## Annotation file

In addition to the reference, we need an annotation file, that tells us, where the genes and proteins are. These are stored in GTF or GFF files. These can ba looked at with a text editor.
[ENSEMBL: homo sapiens](http://www.ensembl.org/Homo_sapiens/Info/Index)

# Trimming data

Often, the FASTQ files still contain adapters, and sometimes low quality bases in the reads (see FASTQC (#TODO)).
Here, we use `trim_galore` to do the adapter trimming. It can automatically run `fastqc` for quality assesment of the FASTQ files.
**For paired end**
Please note that for paired end data, the sequences have to be next to each other in the directory. Also, the `--paired` parameter has to be set.

```r
### Where are the files:
path_to_fq = "/home/raltwas/projekt/data/reads/""
path_to_trimgalore = "/home/raltwas/programme/TrimGalore-0.6.6/trim_galore"

### read all "fastq.gz" files into a list
files_reads <- list.files(path = file.path(path_to_fq),
                          pattern = "*.fastq.gz")

### build the system call
trim_galore_call = paste(path_to_trimgalore,
                         "--fastqc",
                         "-j", "4",
#                        "--paired",
                         "--output_dir", paste0(path_to_fq, "trimmed"),
                         paste(path_to_fq, files_reads, 
                               collapse = " ",
                               sep = ""))

system(trim_galore_call)
```

The result is a copy of all files with the `_trimmed.fq.gz` suffix. These are the trimmed reads.

# Mapping the data

In order to map the FASTQ files to the reference genome, we use `STAR`.

## Indexing

First, `STAR` needs an index of the reference genome. This needs to be done only once for every reference genome. It needs the FASTA of the reference genome, the directory where the index is to be stored, and the GTF file (optional, only for splice junctions. But it's still nice)

```{bash
genome_dir="/data/local/raltwas/index/genomes/Mus_musculus.GRCm38.dna_sm.primary_assembly_cleaned_STARindex"
path_to_reference=""/data/local/raltwas/index/genomes/Mus_musculus.GRCm38.dna_sm.primary_assembly_cleaned.fa"
path_to_gtf="/data/local/raltwas/index/gtf/Mus_musculus.GRCm38.96.gtf"

### create index
STAR \
    --runMode genomeGenerate \
    --genomeDir $genome_dir \
    --genomeFastaFiles $path_to_reference \
    --sjdbGTFfile $path_to_gtf \
    --runThreadN 10 
```

## Mapping

And now the mapping. If the data is single end, you can use the command below. If it is paired end, you need to switch out the `--readFilesIn` line.
***You nee to run this inside the directory of the FASTQ files!**

```{bash
genome_dir="/data/local/raltwas/index/genomes/Mus_musculus.GRCm38.dna_sm.primary_assembly_cleaned_STARindex"
path_to_reference="/data/local/raltwas/index/genomes/Mus_musculus.GRCm38.dna_sm.primary_assembly_cleaned.fa"
path_to_gtf="/data/local/raltwas/index/gtf/Mus_musculus.GRCm38.96.gtf"

### Single End
for i in *_trimmed.fq.gz; do
    STAR \
        --runMode alignReads \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir $genome_dir \
        --sjdbGTFfile path_to_gtf \
        --runThreadN 10 \
        --outFileNamePrefix mapped_reads/${i%_trimmed.fq.gz} \
        --readFilesIn $i
        
### Paired end
   ### --readFilesIn $i ${i%_R1.fastq.gz}_R2.fastq.gz \
done
```

The result is a new directory called `mapped_reads`. Here you find the BAM files. They contain the alignments of the reads and the reference genome.

# DEG analysis

## getting read counts

Here, we read the mapped bam files and collect the counts on all the reads:

```{r
### This reads the differential expression
library("edgeR")
library("Rsubread")
library("rtracklayer")

### read all "bam" files into a list
bam_files <- list.files(path = file.path(mapped_reads),
                          pattern = "*.bam")
path_to_gtf <- "/fast/projects/cubit/current/static_data/annotation/ENSEMBL/75/GRCh37/Homo_sapiens.GRCh37.75.gtf"
path_to_gtf <- "~/projekte/indexes/Homo_sapiens.GRCh38.105.gtf"

### loading Annotation
gtf <- rtracklayer::import(path_to_gtf)
gtf <- as.data.frame(gtf)

### limit to "genes"
gtf <- subset(gtf,
              type == "gene")

### we save the annotation for later
anno <- gtf[, c("gene_id", 
              "gene_name", 
              "gene_biotype")]
save(anno,
     file = "annotation.RData")

### reading of the actual counts
### please note the "isPairedEnd" argument. Adjust accordingly
read_counts <- featureCounts(path_to_bam,
                             annot.ext=gtf,
                             useMetaFeatures=FALSE,
                             allowMultiOverlap=TRUE,
                             tmpDir = tempdir(),
                             isPairedEnd=TRUE,
                             nthreads = 16)

```

## DESeq2

An important part of differential expression analysis is the *sample_table*. It says, which samples are control vs. tumor/treatment. It also writes highlights the replicates.

| sampleName | treatment |
| --- | --- |
| sample_A1 | control |
| sample_A2 | control |
| sample_A3 | control |
| sample_B1 | Knock-down |
| sample_B2 | Knock-down |
| sample_B3 | Knock-down |

```{r volcanoplot}
library("DESeq2")

### !!The order of samples must be the same as in 'feature_counts'!!
sample_table <- data.frame(sampleName = factor(col,
                           treatment =  factor(c("control", "control", "control",
                                                 "knockdown", "knockdown", "knockdown")))

### set the correct control
sample_table[, "treatment"] <- relevel(sample_table[, "treatment"], "control)

### create DESeq data set
dds <- DESeqDataSetFromMatrix(countData = expr_munsch,
                        colData = sample_table,
                        design= ~ treatment)

### do DESeq analysis
dds <- DESeq(dds)

# resultsNames(dds) # lists the coefficients
res <- results(dds, name="treatment_knockdown_vs_control")

### sort results by adjusted p-value
res = res[order(res$padj),]

### get all the significant results
table_sig <- subset(res,
                    padj < 0.05)

table_up <- subset(table_sig,
                   log2FoldChange > 0)
table_down <- subset(table_sig,
                     log2FoldChange < 0)
```

```{r volcanoplot}
 library("ggrastr")    # rasterises the image. Good for images with many dots
library("ggrepel")


color <- matrix("insignificant",
                nrow = nrow(res))
color[res[, "log2FoldChange"] > 0.5] <- "upregulated"
color[res[, "log2FoldChange"] < -0.5] <- "downregulated"

### setting the transparency
alpha <- matrix("insignificant",
                nrow = nrow(res))
alpha[res[, "padj"] < 0.05] <- "significant"

### get the IDs of the 'n' most significant genes
label_ids <- head(order(res[, "padj"]),
                  n = 5)
label <- vector("character", length = nrow(res))
label[label_ids] <- rownames(res)[label_ids]

### pvalue transformation and outlier marking
pval <- log10(res[, "padj"])
shape <- ifelse(pval > 10, "triangle", "circle")
pval[pval > 10] <- 10

### logFC outlier marking
logFC <- res[, "log2FoldChange"]
shape[(abs(logFC) > 4)] <- "triangle"
logFC[logFC >  4] <-  4
logFC[logFC < -4] <- 4

volcano_plot <- ggplot(data.frame(log2FC = logFC,
                                  adj_pValue = pval,
                                  label = label,
                                  shape = shape,
                                  alpha = alpha,
                                  color = color),
                       aes(x = log2FC,
                           y = adj_pValue,
                           label = label,
                           shape = shape,
                           alpha = alpha,
                           color = color)) +

                  scale_color_manual(values = c("upregulated" = "red",
                                                "downregulated" = "blue",
                                                "insignificant" = "grey")) +
                  scale_alpha_manual(guide='none', 
                                     values = c("insignificant"   = 0.3,
                                                "significant" = 1)) +
                  geom_point_rast(size = 2) + # rasterises the image. Good for images with many dots
                  geom_point(size = 2) +
                  xlab(label = "log2FC") +
                  ylab(label = "-log10(adjusted p-value)") +
                  geom_text_repel(size = 5,
                                  point.padding = 0.1,
                                  box.padding = 0.75,
                                  segment.colour = "blue",
                                  colour = "black") +
                  theme(legend.position = "bottom",
                        legend.direction = "horizontal",
                        legend.background = element_rect(color = "steelblue", linetype = "solid"))

print(volcano_plot)
```

