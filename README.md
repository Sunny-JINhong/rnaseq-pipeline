# RNA-seq Pipeline
基于 HISAT2 + StringTie 的 RNA-seq 定量流程（Shell 脚本），自动完成原始测序数据的比对、定量，并生成表达矩阵（count / FPKM / TPM）。

## 流程概述

```
先用HISAT2将RNA-seq数据比对至基因组上，用samtools排序并建索引，用StringTie定量，接着用prepDE.py脚本生成count矩阵，最后提取FPKM/TPM矩阵。
```

单样本处理失败时自动跳过并记录，不中断后续样本；全部完成后统计成功/失败样本数。

## 依赖工具
HISAT2、samtools、StringTie

可通过conda安装：
```
conda install -c bioconda -c conda-forge hisat2 samtools stringtie
```

## 数据准备
需先准备基因组(genome.fasta)及其HISAT2索引和基因组注释文件(genome.gff3)。

## 使用方法

```
bash rnaseq-pipeline.sh <sample_list.txt> <output_dir>
```

## 样本列表格式

制表符(Tab)分隔的三列文本文件，每行一个样本：
```
sample01    /path/to/sample01_R1.fastq.gz    /path/to/sample01_R2.fastq.gz
sample02    /path/to/sample02_R1.fastq.gz    /path/to/sample02_R2.fastq.gz
```

第一列：样本名；第二列：R1 端数据路径；第三列：R2 端数据路径。

## 输出结果
bam文件夹：各样本比对排序后的bam文件及索引；
stringtie文件夹：各样本定量的GTF文件及基因表达表tab文件；
gene_count_matrix.csv：基因 count 矩阵(prepDE.py 生成)；
transcript_count_matrix.csv:转录本 count 矩阵；
FPKM_matrix.csv：基因 FPKM 矩阵；
TPM_matrix.csv：基因 TPM 矩阵；
success_samples.txt/fail_samples.txt：成功/失败样本清单；
run_rnaseq.log：全流程运行日志。

## 注意事项
HISAT2 索引前缀需与脚本中的 `genome.fasta` 一致；
样本列表必须为制表符分隔，列数不足会导致比对报错，样本会记录为失败。
