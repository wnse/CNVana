import logging
import pandas as pd
import os
import time
import argparse

def write_log(logfile, df_dict):
	df = pd.DataFrame.from_dict(df_dict,orient='index')
	df.to_csv(logfile, sep='\t')

def run_ana(inpath, outpath, parralel, threads):
	# inpath = os.path.join(app.config['RAW_PATH'], batch)
	# outpath = os.path.join(app.config["OUT_PATH"], outdir)
	
	if not os.path.exists(outpath):
		os.makedirs(outpath)
	elif os.path.isfile(outpath):
		return 0

	logfile = os.path.join(outpath, 'log')
	status = {}
	sns = []

	for f in os.listdir(inpath):
		suffix = os.path.splitext(f)[1]
		if suffix in ['.fq','.gz','.fastq']:
			file = os.path.join(inpath, f)
			sns.append(file)
			status[file] = 'waiting'
	write_log(logfile, status)
	print(status)
	for s in sns:
		status[s] = 'runing'
		write_log(logfile, status)
		print(status)
		time.sleep(10)

if __name__ == "__main__":
	bin_dir = os.path.split(os.path.realpath(__file__))[0]
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', required=True, help='SE fastq files dir path')
	parser.add_argument('-o', '--outdir', required=True, help='merge bed file and cov file for SCana input')
	parser.add_argument('-t', '--threads', default=2, type=int, help='threads for mapping')
	parser.add_argument('-p', '--parralel', default=2, type=int, help='runing samples simultaneously')
	args = parser.parse_args()


	run_ana(args.input, args.outdir, args.parralel, args.threads)





