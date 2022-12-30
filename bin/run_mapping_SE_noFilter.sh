#!/bin/bash

help()
{
	cat <<- EOF
	Desc: Mapping reads of a genus to UCGs.
	Usage: bash map2UCGs.sh -h
	Author: zhuqh
	Parameters:
		
		-h --help       Print the help info and exit.

		-v --version    Print the version info.

		-i --infastq    The input fastq.

		-o --out		The output files prefix.

		-t --threads    The threads for processing.

		-d --db			The human genome database index.

		-w --win		The human genome windows bed file.

	License: none
	
EOF
	exit 0
}

version_info()
{
	cat <<- EOF

		mapping human genome.sh 2022-6-22

	EOF
		exit 0
}


while true;do
	case "$1" in
		-i|--infq)
			infq=$2;shift 2;;
		-o|--out)
			out=$2;shift 2;;
		-t|--threads)
			threads=$2;shift 2;;
		-d|--db)
			db=$2;shift 2;;
		-w|--win)
			win=$2;shift 2;;
		-h|--help)
			help;;
		-v|--version)
			version_info;;
		*)
			break;
	esac
done

# fastp -w -i /root/PGT_rawdata/HZ20220610/202206082018_210901010_A_PGTA_0608_1_L01.fq.gz -o 202206082018_210901010_A_PGTA_0608_1_L01.clean.fq.gz
# fastp -D -w ${threads} -i ${infq} -o ${out}.clean.fq.gz -j ${out}.fastp.json -h ${out}.fastp.html
# bowtie2 -p 8 -x /root/yangk/db/hg38 -U 202206082018_210901010_A_PGTA_0608_1_L01.clean.fq.gz -S 202206082018_210901010_A_PGTA_0608_1_L01.clean.sam 2>202206082018_210901010_A_PGTA_0608_1_L01.clean.bowtie.log &
bowtie2 -p ${threads} -x ${db} -U ${infq} -S ${out}.sam 2>${out}.bowtie.log
# bowtie2 -p ${threads} -x ${db} -U ${out}.clean.fq.gz -S ${out}.sam 2>${out}.bowtie.log
# samtools view -F 4 -h -q 0 202206082018_210901010_A_PGTA_0608_1_L01.clean.sam |samtools sort - -@ 8 -o 202206082018_210901010_A_PGTA_0608_1_L01.clean.sort.bam
samtools view -F 4 -h -q 2 ${out}.sam | samtools sort - -@ ${threads} -o ${out}.sort.bam
samtools flagstat ${out}.sort.bam >${out}.sort.bam.flagstat
# picard MarkDuplicatesWithMateCigar --REMOVE_DUPLICATES --QUIET -I 202206082018_210901010_A_PGTA_0608_1_L01.clean.sort.bam -O 202206082018_210901010_A_PGTA_0608_1_L01.clean.sort.markdup.bam -M markdup_metrics.txt
picard MarkDuplicatesWithMateCigar --REMOVE_DUPLICATES -I ${out}.sort.bam -O ${out}.sort.markdup.bam -M ${out}.markdup_metrics.txt
samtools index -@ ${threads} ${out}.sort.markdup.bam
samtools idxstats ${out}.sort.markdup.bam >${out}.markdup.inxstats.txt
# bedtools coverage -F 0.5 -a /root/yangk/hg38_max_win.bed -b 202206082018_210901010_A_PGTA_0608_1_L01.clean.sort.markdup.bam
bedtools coverage -F 0.5 -a ${win} -b ${out}.sort.markdup.bam > ${out}.coverage.bed
rm -rf ${out}.sam ${out}.sort.bam ${out}.sort.markdup.bam ${out}.sort.markdup.bam.bai ${out}.clean.fq.gz
