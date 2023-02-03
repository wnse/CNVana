
import os
import re
import json
import argparse
import pandas as pd
from datetime import datetime
from make_report_1 import make_report_1
from make_report_2 import get_config, make_report_2
from plot_chr_scatter import get_plot
from result2exp import dict2ext



def make_PGS_report(df_family, df_sample, inputdir, outdir, config_dir):
	total_sample = df_sample.shape[0]
	report_sample = 0
	temp_types = []
	if '模板分类' in df_family.columns:
		df_family['模板分类'] = df_family['模板分类'].astype(str)
		temp_types = df_family['模板分类'].unique()

	df_sample_info_res = pd.DataFrame()
	if '检测结果' in df_sample.columns:
		res_dict = df_sample.set_index('样本编号')['检测结果'].to_dict()
		df_sample_info_res = dict2ext(res_dict)
		df_sample_info_res = pd.DataFrame.from_dict(df_sample_info_res,orient='index')

		df_sample['结果解释'] = df_sample['样本编号'].map(df_sample_info_res['解释'].to_dict())
		df_sample['是否推荐'] = df_sample['样本编号'].map(df_sample_info_res['推荐'].to_dict())
		# df_sample['结果解释'] = df_sample['结果解释'].replace(';',';\n')

	png_dir = os.path.join(outdir, 'png')
	if not os.path.isdir(png_dir):
		os.makedirs(png_dir)

	sample_csv = {}
	for f in os.listdir(inputdir):
		if os.path.splitext(f)[1] == '.csv':
			if re.match('[A|B]_PGTA_(.+).fq',f):
				sn = re.match('[A|B]_PGTA_(.+).fq',f).group(1)
				sample_csv[sn] = os.path.join(inputdir, f)

	for temp_type in temp_types:
		df_family_tmp = df_family[df_family['模板分类']==temp_type].fillna('').copy()
		df_family_tmp['index'] = df_family_tmp['家系编号']
		df_sample_tmp = df_sample[df_sample['家系编号'].isin(df_family_tmp['家系编号'].to_list())].fillna('').copy()
		df_sample_tmp['index'] = df_sample[['家系编号','样本编号']].apply(lambda x: ':'.join(x.to_list()), axis=1)


		if temp_type != 'nan':
			temp_num, fig_num = str(temp_type).split('.')
			for idx in df_sample_tmp['index']:
				f_idx, sn = idx.split(":")
				sex_chr = ''
				if '性染色体' in df_sample_tmp.columns:
					sex_chr = df_sample_tmp.loc[df_sample_tmp[df_sample_tmp['index']==idx].index[0],'性染色体']
				report_sample += 1
				if sn in sample_csv.keys():
					# print(sample_csv[sn], fig_num, png_dir)
					png_file = os.path.join(png_dir, f'{sn}.{fig_num}.png')
					cmd = f"python {os.path.join(bin_dir, 'plot_chr_scatter.py')} -i {sample_csv[sn]} -o {png_file} -t {fig_num} -s {sex_chr}"
					try:
						cmd_out = os.system(cmd)
					except Exception as e:
						cmd_out = e
					with open(logfile, 'a') as h:
						print(f"{sn}\tmakePNG\t{cmd_out}", file=h, flush=True)	
				else:
					with open(logfile, 'a') as h:
						print(f"{sn}\tmakePNG ERROR\t无csv", file=h, flush=True)			
			
			f_config = os.path.join(config_dir, 'template_config.xlsx')
			
			temp_docx = os.path.join(config_dir, f'template_{temp_num}.docx')
			if '送检区域' in df_family_tmp.columns:
				for region in df_family_tmp['送检区域'].unique():
					df_family_tmp_region = df_family_tmp[df_family_tmp['送检区域']==region]
					df_sample_tmp_region = df_sample_tmp[df_sample_tmp['家系编号'].isin(df_family_tmp_region['家系编号'].to_list())]
					outdir_region = os.path.join(outdir, region)
					if not os.path.isdir(outdir_region):
						os.makedirs(outdir_region)
					if temp_num == '1':
						dict_family = df_family_tmp.set_index('index').to_dict(orient='index')
						dict_sample = df_sample_tmp.groupby('家系编号').apply(lambda x: x.set_index('样本编号').to_dict(orient='index'))
						df_config = pd.read_excel(f_config, f'{temp_num}', index_col=0)
						dict_config = df_config.to_dict(orient='index')
						make_report_1(temp_docx, dict_family, dict_sample, dict_config, outdir=outdir_region, png_dir=png_dir, png_name=int(fig_num))

					if temp_num in ['2','3']:
						df_config = get_config(f_config, f'{temp_num}')
						png_suffix = f'.{str(fig_num)}.png'
						make_report_2(df_family_tmp.set_index('index'), df_sample_tmp.set_index('index'), df_config, temp_docx, outdir_region, png_dir, png_suffix)
			else:
				if temp_num == '1':
					dict_family = df_family_tmp.set_index('index').to_dict(orient='index')
					dict_sample = df_sample_tmp.groupby('家系编号').apply(lambda x: x.set_index('样本编号').to_dict(orient='index'))
					df_config = pd.read_excel(f_config, f'{temp_num}', index_col=0)
					dict_config = df_config.to_dict(orient='index')
					make_report_1(temp_docx, dict_family, dict_sample, dict_config, outdir=outdir, png_dir=png_dir, png_name=int(fig_num))

				if temp_num in ['2','3']:
					df_config = get_config(f_config, f'{temp_num}')
					png_suffix = f'.{str(fig_num)}.png'
					make_report_2(df_family_tmp.set_index('index'), df_sample_tmp.set_index('index'), df_config, temp_docx, outdir, png_dir, png_suffix)

			# dt = datetime.strftime(datetime.now(),format='%Y-%m-%d %H:%M:%S')
			# with open(logfile, 'a') as h:
			# 	print(f"{temp_type}\tDONE\t{dt}", file=h, flush=True)
		else:
			for i in df_family[df_family['模板分类'] == 'nan']['家系编号'].to_list():
				with open(logfile, 'a') as h:
					print(f"{i}\tERROR\t模板错误", file=h, flush=True)
	dt = datetime.strftime(datetime.now(),format='%Y-%m-%d %H:%M:%S')
	with open(logfile, 'a') as h:
		print(f"{report_sample}/{total_sample}\tDONE\t{dt}", file=h, flush=True)


def check_sample_info(config_dir, batch_dir, logfile=None):
	df_family_info = pd.DataFrame()
	with open(os.path.join(config_dir,'template_type.json'),'rt') as h:
		template_all = json.load(h)
	file = os.path.join(batch_dir, 'sample_info.xlsx')
	if os.path.isfile(file):
		df_family_info = pd.read_excel(file, '家系')
		if '模板' in df_family_info.columns:
			df_family_info['模板分类'] = df_family_info['模板'].map(template_all)
		df_sample_info = pd.read_excel(file, '样本')
	else:
		if logfile:
			with open(logfile, 'a') as h:
				print(f"文件错误\tERROR\t没有文件【{file}】", file=h, flush=True)

	return df_family_info, df_sample_info

if __name__ == '__main__':
	bin_dir = os.path.split(os.path.realpath(__file__))[0]
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--indir', required=True, help='input dir for excel and csv')
	parser.add_argument('-c', '--config', default=bin_dir, help='config dir')
	parser.add_argument('-o', '--outdir', default=None, help='out dir for report')
	args = parser.parse_args()

	inputdir = args.indir
	config_dir = args.config
	outdir = args.outdir
	
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
	logfile = os.path.join(outdir,'log')

	df_family_info, df_sample_info = check_sample_info(config_dir, inputdir, logfile)
	make_PGS_report(df_family_info, df_sample_info, inputdir, outdir, config_dir)

