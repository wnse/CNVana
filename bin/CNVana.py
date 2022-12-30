# -*- coding: utf-8 -*-
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
import sys
import shutil
import argparse
import pandas as pd

from mapping import fastp_res, sta_coverage, format_input_for_SCana


# %%
def run_mapping(bindir, fq, threads, db, bed_file, outdir, sn):
    tmp_dir = os.path.join(outdir, sn)
    cmd = f"sh {os.path.join(bindir, 'run_mapping_SE_noFilter.sh')} -i {fq} -o {tmp_dir} -t {threads} -d {db} -w {bed_file}"
    print(cmd)
    cov_file = os.path.join(outdir, sn+".coverage.bed")
    if os.path.exists(cov_file):
        return f'{cov_file} exists！skip mapping'
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    cmd_out = os.system(cmd)
    return cmd_out



# %%
def run_mapping_res(gc_file, json_file, cov_file, SCana_input_file):
    mapping_out = pd.DataFrame()
    fastp_out = pd.DataFrame()
    if os.path.isfile(json_file):
        fastp_out = fastp_res(json_file)
    cover_sta = pd.DataFrame()
    if os.path.isfile(cov_file):
        outdir = os.path.split(SCana_input_file)[0]
        cover_sta = sta_coverage(cov_file)
        format_input_for_SCana(gc_file, cov_file, SCana_input_file)

        mapping_out = pd.concat([fastp_out, cover_sta])
        mapping_out.to_csv(os.path.join(outdir,'sta.csv'),header=None) 
    return mapping_out


# %%
def run_SCana(bindir, SCana_input_file, BG_file, outdir, sn):
    if not os.path.isfile(SCana_input_file):
        return f'{SCana_input_file} not exists！skip SCana'
#     bg = os.path.splitext(os.path.split(BG_file)[1])[0]
#     tmp_dir = os.path.join(outdir, bg)
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    cmd = f"Rscript {os.path.join(bindir, 'SC_Analysis_DNAcopy.r')} {outdir} {SCana_input_file} {BG_file} {sn}"
    print(cmd)
    cmd_out = os.system(cmd)
    return cmd_out



# %%
def run_SCana_res(SCana_summary_file, outdir, sn, tag):
    df = pd.read_csv(SCana_summary_file)[['SD','mean','sexchr']].loc[0]
    df['cv'] = df['SD']/df['mean']
    df.index = df.index + "_" + str(tag)
    try:
        out_file = os.path.join(outdir, 'sta.csv')
        df_sta = pd.read_csv(out_file,header=None, index_col=0)
        pd.concat([df_sta[1], df]).to_csv(out_file, header=None)
    except Exception as e:
        df.to_csv(out_file, header=None)


# %%
if __name__ == "__main__":
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='SE fastq file')
    parser.add_argument('-n', '--name', required=True, help='sample name')
    parser.add_argument('-o', '--outdir', required=True, help='merge bed file and cov file for SCana input')
    parser.add_argument('-d', '--db', required=True, help='human genome index for bowtie')
    parser.add_argument('-t', '--threads', default=2, help='threads for mapping')
    parser.add_argument('-b', '--bed_file', default=os.path.join(bin_dir, 'bedfile','hg19_genemind.bed'), help='bed file for SCana including GC Map etc')
    parser.add_argument('-g', '--gc_file', default=os.path.join(bin_dir, 'bedfile','SC_hg19_bin150k.bed'), help='bed file for SCana including GC Map etc')
    parser.add_argument('-bg', '--bg_file', nargs='+', default=[os.path.join(bin_dir,'baseline', 'G1129raw.txt')], help='baseline file for correction')
    args = parser.parse_args()   
    
    
    cmd_out = run_mapping(bin_dir, args.input, args.threads, args.db, args.bed_file, args.outdir, args.name)
    print(f'run_mapping: {cmd_out}')
    json_file = os.path.join(args.outdir, args.name + ".fastp.json")
    cov_file = os.path.join(args.outdir, args.name + ".coverage.bed")
    SCana_input_file = os.path.join(args.outdir, args.name + ".gm.cov.bed")
    run_mapping_res(args.gc_file, json_file, cov_file, SCana_input_file)
    
    for bg_file  in args.bg_file:
        bg = os.path.splitext(os.path.split(bg_file)[1])[0]
        tmp_dir = os.path.join(args.outdir, bg)
        SCana_summary_file = os.path.join(tmp_dir, args.name + "_summary.csv")
        cmd_out = run_SCana(bin_dir, SCana_input_file, bg_file, tmp_dir, args.name)
        run_SCana_res(SCana_summary_file, args.outdir, args.name, bg)
        print(f'run_SCana: {cmd_out}')
    os.remove(SCana_input_file)

