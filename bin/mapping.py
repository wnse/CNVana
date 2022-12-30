# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import os
import json
# import logging
import argparse
import pandas as pd

# import matplotlib.pyplot as plt
# import seaborn as sns

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# ## step
# ```
# 1. mapping
# 2. generate_bedfile
# 3. sc_analysis
# 4. generate_single_res
#     1. sta
#     2. cnv res
#     3. cov & plot
# 5. generate_batch_res
# ```

# %% [markdown]
# ## mapping

# %% [markdown] tags=[]
#
# ```
# fastp -D -w ${threads} -i ${infq} -o ${out}.clean.fq.gz -j ${out}.fastp.json -h ${out}.fastp.html
# bowtie2 -p ${threads} -x ${db} -U ${out}.clean.fq.gz -S ${out}.sam 2>${out}.bowtie.log
# samtools view -F 4 -h -q 2 ${out}.sam | samtools sort - -@ ${threads} -o ${out}.sort.bam
# samtools flagstat ${out}.sort.bam >${out}.sort.bam.flagstat
# picard MarkDuplicatesWithMateCigar --REMOVE_DUPLICATES -I ${out}.sort.bam -O ${out}.sort.markdup.bam -M ${out}.markdup_metrics.txt
# bedtools coverage -F 0.5 -a ${win} -b ${out}.sort.markdup.bam > ${out}.coverage.bed
#
# ```

# %% [markdown]
# ```
# docker run --rm -v /Users/yangkai/work/GermountX/PGT/20221219_pipline/data:/data/ cnv_ana_build_test:1.2 sh /data/out1/run.sh
#
# python /data/mapping.py -o /data/out/PGTA_GLP-11_gm.cov.bed -c /data/out/PGTA_GLP-11.coverage.bed -j /data/out/PGTA_GLP-11.fastp.json
# Rscript /data/SC_Analysis_DNAcopy.r /data/out PGTA_GLP-11_gm.cov.bed /data/db/G1129raw.txt PGTA_GLP-11_G1129
# Rscript /data/SC_Analysis_DNAcopy.r /data/out PGTA_GLP-11_gm.cov.bed /data/db/BG_all2.txt PGTA_GLP-11_BG_all2
#
# ```

# %%
def fastp_res(json_file):
    base_sta = {}
    tmp = None
    if os.path.isfile(json_file):
        with open(json_file, 'rt') as h:
            tmp = json.load(h)
    base_sta = pd.DataFrame(tmp['summary'])[['before_filtering','after_filtering']].unstack()
    base_sta.index = [i[0].replace('before_filtering','raw').replace('after_filtering','clean')
                      + "_"+i[1] 
                      for i in base_sta.index.to_list()]
    base_sta['duplication'] = tmp['duplication']['rate']
    base_sta = pd.concat([base_sta,pd.Series(tmp['filtering_result'])])
    return base_sta
    


# %%
def sta_coverage(cov_file):
    df_cov = pd.read_csv(cov_file,sep='\t',header=None)
    df_cov.columns = ['chr','start','end','reads','bases','length','cover']
    df_cov['cover'].describe()
    
    chr_list = ['chr'+str(i) for i in range(1,23)]
    cov_sta = df_cov.loc[df_cov['chr'].isin(chr_list),'cover'].describe().round(3)
    
#     cov_sta = [str(i)+':'+ str(v) for i,v in cov_sta.to_dict().items()]
#     sns.boxplot(x='chr',y='cover',data=df_cov, color='white', flierprops={"marker": "x"})
#     plt.annotate(text='|'.join(cov_sta[:2]), xy=(0,0.18), xytext=(0,0.18))
#     plt.annotate(text='|'.join(cov_sta[2:]), xy=(0,0.17), xytext=(0,0.17))
#     plt.ylim(0,0.2)
#     plt.xticks(rotation=90)
#     plt.hlines(y=df_cov['cover'].mean(), xmin=0, xmax=23, color='red', alpha=0.5)
#     plt.hlines(y=0.04, xmin=0, xmax=23, color='grey', alpha=0.5)
#     plt.show()

    cov_sta = cov_sta.drop('count')
    cov_sta.index = 'cover_' + cov_sta.index
    return cov_sta


# %%
def format_input_for_SCana(bed_file, cov_file, out_file):
    df_cov = pd.read_csv(cov_file, sep='\t', header=None, index_col=[0,1])
    df_bed = pd.read_csv(bed_file, sep='\t', header=None, index_col=[0,1])
    df_out = pd.concat([df_bed, df_cov[3]], axis=1)
    df_out.to_csv(out_file, header=None)


# %%
if __name__ == "__main__":
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--out_file', required=True, help='merge bed file and cov file for SCana input')
    parser.add_argument('-c', '--cov_file', default=None, help='bedtools coverage sta file')
    parser.add_argument('-j', '--json_file', default=None, help='fastp json file')
    parser.add_argument('-b', '--bed_file', default=os.path.join(bin_dir, 'bedfile','GM_hg19_bin150k.bed'), help='bed file for SCana including GC Map etc')
    args = parser.parse_args()          

    json_file = args.json_file
    cov_file = args.cov_file
    bed_file = args.bed_file
    SCana_input_file = args.out_file
    outdir = os.path.split(SCana_input_file)[0]

    fastp_out = pd.DataFrame()
    if json_file:
        fastp_out = fastp_res(json_file)

    cover_sta = pd.DataFrame()
    if cov_file:
        cover_sta = sta_coverage(cov_file)
        format_input_for_SCana(bed_file, cov_file, SCana_input_file)

        mapping_out = pd.concat([fastp_out, cover_sta])
        mapping_out.to_csv(os.path.join(outdir,'sta.csv'),header=None)

# %%
# sn = 'PGTA_GLP-11'
    
# json_file = '../data/out/PGTA_GLP-11.fastp.json'
# cov_file = '../data/out/PGTA_GLP-11.coverage.bed'
# bed_file = 'bedfile/GM_hg19_bin150k.bed'
# SCana_input_file = '../data/out1/PGTA_GLP-11_gm.cov.bed'

# fastp_out = fastp_res(json_file)
# cover_sta = sta_coverage(cov_file)
# mapping_out = pd.concat([fastp_out, cover_sta])

# mapping_out.to_csv(sn+'_base_sta.csv',header=None)
# format_input_for_SCana(bed_file, cov_file, SCana_input_file)

# %%

# %%
