#!/bin/bash

# Usage：bash run_rnaseq.sh sample_list.txt output_dir

if [ "$#" -ne 2 ]; then
    echo "Usage: bash run_rnaseq.sh <sample_list.txt> <output_dir>"
    exit 1
fi

SAMPLE_LIST=$1 
OUTPUT_DIR=$2

mkdir -p ${OUTPUT_DIR}/bam ${OUTPUT_DIR}/stringtie ${OUTPUT_DIR}/logs

RUN_LOG="${OUTPUT_DIR}/run_rnaseq.log"

cd ${OUTPUT_DIR}

exec > >(tee -a ${RUN_LOG})
exec 2>&1

echo "===== [$(date)] 开始运行 ====="
echo "样品列表: ${SAMPLE_LIST}"
echo "输出目录: ${OUTPUT_DIR}"
echo ""

SUCCESS_LIST="${OUTPUT_DIR}/success_samples.txt"
FAIL_LIST="${OUTPUT_DIR}/fail_samples.txt"

> "${SUCCESS_LIST}"
> "${FAIL_LIST}"

while IFS=$'\t' read -r sample r1 r2 || [ -n "${sample}" ]; do

	echo "===== [$(date)] 处理样品: ${sample} ====="

	# ---------- 步骤1：HISAT2 比对 ----------

	hisat2 -p 8 --dta -x genome.fasta -1 ${r1} -2 ${r2} -S ${OUTPUT_DIR}/${sample}_vs_genome.sam 2> ${OUTPUT_DIR}/logs/${sample}_vs_genome.log

	if [ $? -ne 0 ] || [ ! -f "${sample}_vs_genome.sam" ] || [ ! -s "${sample}_vs_genome.sam" ]; then
		echo "[$(date)] ${sample} HISAT2 比对失败，跳过"
		rm -f ${sample}_vs_genome.sam
		echo "${sample}" >> ${FAIL_LIST}
		echo ""
		continue
	fi

	samtools sort -@ 4 -o ${OUTPUT_DIR}/bam/${sample}_vs_genome.bam ${OUTPUT_DIR}/${sample}_vs_genome.sam

	if [ $? -ne 0 ] || [ ! -f "bam/${sample}_vs_genome.bam" ] || [ ! -s "bam/${sample}_vs_genome.bam" ]; then
		echo "[$(date)] ${sample} samtools sort 失败，跳过"
		rm -f ${sample}_vs_genome.sam bam/${sample}_vs_genome.bam
		echo "${sample}" >> ${FAIL_LIST}
		echo ""
		continue
	fi

	samtools index ${OUTPUT_DIR}/bam/${sample}_vs_genome.bam

	rm -f ${OUTPUT_DIR}/${sample}_vs_genome.sam

	# ---------- 步骤2：StringTie 定量 ----------

	stringtie -e -B -p 8 -G genome.gff3 -o ${OUTPUT_DIR}/stringtie/${sample}_vs_genome.gtf -A ${OUTPUT_DIR}/stringtie/${sample}_vs_genome.tab ${OUTPUT_DIR}/bam/${sample}_vs_genome.bam

	if [ $? -ne 0 ] || [ ! -f "stringtie/${sample}_vs_genome.gtf" ] || [ ! -s "stringtie/${sample}_vs_genome.gtf" ]; then
		echo "[$(date)] ${sample} StringTie 定量失败，跳过"
		echo ${sample} >> ${FAIL_LIST}
		echo ""
		continue
	fi

	echo "[$(date)] ${sample} 处理成功"
	echo ${sample} >> ${SUCCESS_LIST}
	echo "===== [$(date)] ${sample} 比对完成 ====="
	echo ""

done < ${SAMPLE_LIST}

echo "===== [$(date)]所有样品比对和定量完成 ====="
echo "成功: $(wc -l < ${SUCCESS_LIST}) 个"
if [ -s ${FAIL_LIST} ]; then
	echo "失败: $(wc -l < ${FAIL_LIST}) 个"
	echo "失败样品:"
	cat ${FAIL_LIST}
fi
echo ""

# ---------- 步骤3：prepDE.py ----------

echo "[$(date)] 生成 mergelist.txt..."

> ${OUTPUT_DIR}/mergelist.txt

while IFS=$'\t' read -r sample r1 r2 || [ -n "${sample}" ]; do
	if grep -q "^${sample}$" ${SUCCESS_LIST}; then
		echo -e "${sample}\t${OUTPUT_DIR}/stringtie/${sample}_vs_genome.gtf" >> ${OUTPUT_DIR}/mergelist.txt
	fi
done < ${SAMPLE_LIST}

if [ ! -s mergelist.txt ]; then
	echo "[$(date)] 没有成功处理的样品，无法运行 prepDE.py"
        exit 1
fi

echo "[$(date)] 运行 prepDE.py 提取表达矩阵..."

prepDE.py -i mergelist.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv 2> prepDE.log

echo "[$(date)] prepDE.py 完成"
echo "基因表达矩阵: ${OUTPUT_DIR}/gene_count_matrix.csv"
echo "转录本表达矩阵: ${OUTPUT_DIR}/transcript_count_matrix.csv"

echo ""

# ---------- 提取 FPKM/TPM 矩阵 ----------
echo "[$(date)] 提取 FPKM/TPM 矩阵..."
if [ ! -s ${SUCCESS_LIST} ]; then
	echo "[$(date)] 没有成功处理的样品，无法提取 FPKM/TPM"
else
	echo "[$(date)] 成功样品数: $(wc -l < ${SUCCESS_LIST})"

	TAB_FILES=""
	while read -r sample; do
        tab_file="stringtie/${sample}_vs_genome.tab"
        if [ ! -f "${tab_file}" ]; then
            echo "[$(date)] ⚠️  文件不存在，跳过: ${tab_file}"
            continue
        fi
        if [ -z "${TAB_FILES}" ]; then
            TAB_FILES="${tab_file}"
        else
            TAB_FILES="${TAB_FILES} ${tab_file}"
        fi
    done < ${SUCCESS_LIST}

    if [ -z "${TAB_FILES}" ]; then
        echo "[$(date)] ❌ 没有可用的 tab 文件"
    else
        # 提取 FPKM 矩阵
        {
            echo -n "Gene_ID"
            while read -r sample; do
                echo -n ",${sample}"
            done < ${SUCCESS_LIST}
            echo ""
            
            paste ${TAB_FILES} | awk 'NR>1 {
                printf "%s", $1;
                for(i=8; i<=NF; i+=9) printf ",%s", $i;
                printf "\n"
            }'
        } > FPKM_matrix.csv

        # 提取 TPM 矩阵
        {
            echo -n "Gene_ID"
            while read -r sample; do
                echo -n ",${sample}"
            done < ${SUCCESS_LIST}
            echo ""
            
            paste ${TAB_FILES} | awk 'NR>1 {
                printf "%s", $1;
                for(i=9; i<=NF; i+=9) printf ",%s", $i;
                printf "\n"
            }'
        } > TPM_matrix.csv

        echo "[$(date)] FPKM/TPM 矩阵提取完成"
        echo "FPKM 矩阵: FPKM_matrix.csv"
        echo "TPM 矩阵: TPM_matrix.csv"
    fi
fi

echo "===== [$(date)] 全部完成！ ====="
