
# 🧬 Variant Calling Workflow using BWA, SAMtools, and BCFtools

This guide outlines the steps to identify genetic variants in *E. coli* using publicly available tools and data. It is meant for educational use, particularly for school students learning about bioinformatics pipelines.

---

## 📌 Workflow Overview


---

## 🛠️ Step 1: Install Required Tools

```bash
sudo apt install bwa samtools bcftools
```

Tools installed:
- `bwa`: For aligning sequencing reads.
- `samtools`: For handling SAM/BAM files.
- `bcftools`: For variant calling.

---

## 📂 Step 2: Download Reference Genome

```bash
mkdir -p /home/data/ref_genome

curl -L -o /home/data/ref_genome/ecoli_rel606.fasta.gz \
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz

gunzip /home/data/ref_genome/ecoli_rel606.fasta.gz
```

---

## 📦 Step 3: Download Sample Read Data

```bash
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
mv sub/ /home/data/trimmed_fastq_small
```

This dataset contains paired-end FASTQ files.

---

## 🧬 Step 4: Index the Reference Genome

```bash
bwa index /home/data/ref_genome/ecoli_rel606.fasta
```

Indexing prepares the reference for alignment.

---

## 🧾 Step 5: Align Reads with BWA

```bash
bwa mem /home/data/ref_genome/ecoli_rel606.fasta \
data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq \
data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq \
> /home/data/SRR2584866.aligned.sam
```

Output: `SRR2584866.aligned.sam`

---

## 🔄 Step 6: Convert SAM to BAM

```bash
samtools view -S -b /home/data/SRR2584866.aligned.sam > /home/data/SRR2584866.aligned.bam
```

BAM files are binary and compressed.

---

## 📚 Step 7: Sort and index the BAM File

```bash
samtools sort -o /home/data/SRR2584866.aligned.sorted.bam /home/data/SRR2584866.aligned.bam
```

```bash
 samtools index SRR2584866.aligned.sorted.bam
```
Sorting and indexing improve processing speed for the next steps.

---

## 📊 Step 8: Check BAM File Stats

```bash
samtools flagstat /home/data/SRR2584866.aligned.sorted.bam
```

This gives summary metrics such as total reads and mapping rates.

---

## 🔍 Step 9: Generate Raw BCF File

```bash
bcftools mpileup -O b -o /home/data/SRR2584866_raw.bcf \
-f /home/data/ref_genome/ecoli_rel606.fasta \
/home/data/SRR2584866.aligned.sorted.bam
```

---

## 🧬 Step 10: Call Variants

```bash
bcftools call --ploidy 1 -m -v -o /home/data/SRR2584866_variants.vcf /home/data/SRR2584866_raw.bcf
```

- `--ploidy 1`: Because bacteria are haploid.
- Output: VCF file listing SNPs/indels.

---

## 👁️ Step 11: View Variants

```bash
less -S /home/data/SRR2584866_variants.vcf
```

Use arrow keys to scroll horizontally if needed.

---

## ✅ Conclusion

This workflow covers a full variant calling pipeline starting from raw reads through to variant identification. It introduces essential tools and formats (FASTQ, FASTA, SAM, BAM, BCF, VCF) in genome analysis.

---
