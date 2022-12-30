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
import json
import argparse
import shutil
import pandas as pd


# %%
import multiprocessing
from datetime import datetime

# %%
from CNVana import run_mapping, run_mapping_res, run_SCana, run_SCana_res


# %%
def CNVana_batch(fastq, db, threads, outdir, name, gc_file, bed_file, bg_files):
    bg_out_file = {}
    try:
        dt = datetime.strftime(datetime.now(),format='%Y-%m-%d %H:%M:%S')
        with open(os.path.join(outdir, 'log'),'a') as h:
            print(f"{name}\trun\t{dt}", file=h, flush=True)
            
        cmd_out = run_mapping(bin_dir, fastq, threads, db, bed_file, os.path.join(outdir,name), name)
        print(f'{name}:run_mapping: {cmd_out}')
        json_file = os.path.join(outdir, name, name + ".fastp.json")
        cov_file = os.path.join(outdir, name, name + ".coverage.bed")
        SCana_input_file = os.path.join(outdir, name, name + ".gm.cov.bed")
        run_mapping_res(gc_file, json_file, cov_file, SCana_input_file)


        for bg_file  in bg_files:
            bg = os.path.splitext(os.path.split(bg_file)[1])[0]
            tmp_dir = os.path.join(outdir, name, bg)
            SCana_summary_file = os.path.join(tmp_dir, name + "_summary.csv")
            cmd_out = run_SCana(bin_dir, SCana_input_file, bg_file, tmp_dir, name)
            run_SCana_res(SCana_summary_file, os.path.join(outdir,name), name, bg)
            print(f'run_SCana: {cmd_out}')
            bg_out_file[bg] = tmp_dir
        os.remove(SCana_input_file)
        dt = datetime.strftime(datetime.now(),format='%Y-%m-%d %H:%M:%S')
        with open(os.path.join(outdir, 'log'),'a') as h:
            print(f"{name}\tDONE\t{dt}", file=h, flush=True)
    except Exception as e:
        print(e)
    return bg_out_file
    
def CNVana_batch_res(samples, stas, bg_out_files, outdir):
    df_sta = pd.DataFrame()

    for i, s in enumerate(samples):
        df_tmp = pd.read_csv(stas[i], header=None, index_col=0)
        df_tmp.columns = [s]
        df_sta = pd.concat([df_sta, df_tmp],axis=1)
        for bg in bg_out_files[i]:
            batch_bg_dir = os.path.join(outdir,bg)
            if not os.path.isdir(batch_bg_dir):
                os.makedirs(batch_bg_dir)
            for f in os.listdir(bg_out_files[i][bg]):
                shutil.move(os.path.join(bg_out_files[i][bg],f), batch_bg_dir)
    df_sta.to_csv(os.path.join(outdir,'batch_sta.csv'))

def run_CNVana_batch(indir, db, threads, outdir, split_pattern='', start=3, end=-1):
    samples = []
    stas = []
    bg_out_files = []

    for f in os.listdir(indir):
        suffix = os.path.splitext(f)[1]
        if suffix in ['.fq','.gz','.fastq']:
            print(f"{name:=^100}")
            fastq = os.path.join(indir, f)
            if split_pattern:
                name = '_'.join(f.split(split_pattern)[start:end])
            else:
                name = suffix = os.path.splitext(f)[0]
            print(f"fastp{name:=^100}")
            try:
                bg_out_file = CNVana_batch(fastq, db, threads, outdir, name, gc_file, bed_file, bg_files)
            except Exception as e:
                bg_out_file = []
                print(f"CNVana_batch error {name}: {e}")
                
            sta_file = os.path.join(outdir, name, 'sta.csv')
            if os.path.isfile(sta_file):
                stas.append(sta_file)
                bg_out_files.append(bg_out_file)
                samples.append(name)
    CNVana_batch_res(samples, stas, bg_out_files, os.path.join(outdir,'CNVana_batch_result'))
            
    

# %%
def run_fp(infq, outdir, name, threads):
    outdict = {}
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfq = os.path.join(outdir, f"{name}.clean.fq.gz")
    outjson = os.path.join(outdir, f"{name}.fastp.json")
    outhtml = os.path.join(outdir, f"{name}.fastp.html")
    cmd = f"fastp -D -i {infq} -o {outfq} -j {outjson} -h {outhtml}"
    print(f"{name}:{cmd}")
    cmd_out = os.system(cmd)
    outdict[name] = outjson
    return outfq



# %%
# def run_test(infq, outdir, name, threads):
#     print(f"{name}:{infq}")
#     time.sleep(2)
#     print(f"{name}:{outdir}")
    


# %%
def run_CNVana_batch_multi(indir, db, threads, outdir, split_pattern='', start=3, end=-1, multi=2):
    samples = []
    infiles = []
    stas = []
    bg_out_files = []
    pool_res = []
    total_sample = 0

    
    for f in os.listdir(indir):
        suffix = os.path.splitext(f)[1]
        if suffix in ['.fq','.gz','.fastq']:
            if split_pattern:
                name = split_pattern.join(f.split(split_pattern)[start:end])
                if not name:
                    name = os.path.splitext(f)[0]
            else:
                name = os.path.splitext(f)[0]
            samples.append(name)
            infiles.append(f)
            dt = datetime.strftime(datetime.now(),format='%Y-%m-%d %H:%M:%S')
            with open(os.path.join(outdir, 'log'),'a') as h:
                print(f"{name}\twait\t{dt}", file=h, flush=True)      

    pool = multiprocessing.Pool(processes=multi)
    for i, name in enumerate(samples):
        f = infiles[i]
        fastq = os.path.join(indir, f)
        try:
            filter_fq = run_fp(fastq, os.path.join(outdir,name), name, threads)
        except Exception as e:
            bg_out_file = []
            print(f"CNVana_batch filter_fq error {name}: {e}")
        try:
            pool_res.append(pool.apply_async(CNVana_batch, args=(filter_fq, db, threads, outdir, name, gc_file, bed_file, bg_files,)))
        except Exception as e:
            bg_out_file = []
            print(f"CNVana_batch multi error {name}: {e}")
        sta_file = os.path.join(outdir, name, 'sta.csv')
        stas.append(sta_file)
            
    pool.close()
    pool.join()
    for i in pool_res:
        print(i.get())             
        bg_out_files.append(i.get())
    CNVana_batch_res(samples, stas, bg_out_files, os.path.join(outdir,'CNVana_batch_result'))



# %%
if __name__ == "__main__":
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='SE fastq files dir path')
#     parser.add_argument('-n', '--name', required=True, help='sample name')
    parser.add_argument('-o', '--outdir', required=True, help='merge bed file and cov file for SCana input')
    parser.add_argument('-d', '--db', required=True, help='human genome index for bowtie')
    parser.add_argument('-t', '--threads', default=2, help='threads for mapping')
    parser.add_argument('-b', '--bed_file', default=os.path.join(bin_dir, 'bedfile','hg19_genemind.bed'), help='bed file for mapping coverage')
    parser.add_argument('-g', '--gc_file', default=os.path.join(bin_dir, 'bedfile','SC_hg19_bin150k.bed'), help='bed file for SCana including GC Map etc')
    parser.add_argument('-bg', '--bg_file', action='append', default=[], help='baseline file for correction')
    parser.add_argument('-p', '--pattern', default='', help='fastq file name split to sample name by pattern')
    parser.add_argument('-m', '--multi', default=1, type=int, help='parralel samples for analysis')

    args = parser.parse_args()
    gc_file, bed_file, bg_files = args.gc_file, args.bed_file, args.bg_file
#     if args.multi:
    run_CNVana_batch_multi(args.input, args.db, args.threads, args.outdir, split_pattern=args.pattern, multi=args.multi)
#     else:
#         run_CNVana_batch(args.input, args.db, args.threads, args.outdir, split_pattern=args.pattern)

    
    

