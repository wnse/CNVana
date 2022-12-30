
import os
import shutil
import subprocess
import multiprocessing
import flask
import json
import pandas as pd
from flask import render_template, session, redirect, url_for, request, jsonify
from cnv_ana import app
from cnv_ana.forms import BasesLookupForm, AnaForm
import multiprocessing


class ChoiceObj(object):
    def __init__(self, name, choices):
        # this is needed so that BaseForm.process will accept the object for the named form,
        # and eventually it will end up in SelectMultipleField.process_data and get assigned
        # to .data
        setattr(self, name, choices)

@app.route("/")
def home():
    return render_template('base.html')


def get_batches():
    input_list = []
    raw_path = app.config['RAW_PATH']
    for d in os.listdir(raw_path):
        if os.path.isdir(os.path.join(raw_path,d)):
            input_list.append(d)
    input_list = sorted(input_list, key=lambda x: os.path.getmtime(os.path.join(raw_path, x)))
    return input_list

def get_baseline():
    baseline_list = []
    base_path = app.config['BASE_PATH']
    for b in os.listdir(base_path):
        if os.path.isfile(os.path.join(base_path,b)):
            baseline_list.append(b.strip('.txt'))
    return baseline_list


def check_info(df):
    global cpu_count
    info = df['信息确认']
    check = {}
    outpath = app.config['RAW_PATH']
    outdir = os.path.join(outpath,info['批次'])
    if os.path.exists(outdir):
        if os.path.isdir(outdir):
            check['批次'] = 'OK'
        else:
            check['批次'] = f"【{info['批次']}】不是一个目录，需要修改批次"


    outpath = app.config['OUT_PATH']
    outdir = os.path.join(outpath,info['输出'])
    if os.path.exists(outdir):
        if os.path.isdir(outdir):
            check['输出'] = f"存在目录【{info['输出']}】，最好修改输出目录。"
        else:
            check['输出'] = f"存在文件【{info['输出']}】，需要修改输出目录。"
    else:
        check['输出'] = 'OK'

    basepath = app.config['BASE_PATH']
    tmp = []
    for b in info['基线'].split(';'):
        bf = os.path.join(basepath, b+'.txt')
        if not os.path.isfile(bf):
            tmp.append(f'【{bf}】不存在，请返回重新选择基线。')
    if tmp:
        check['基线'] = ';'.join(tmp)
    else:
        check['基线'] = 'OK'

    
    if info['并行'] > cpu_count:
        check['并行'] = f"CPU总核心数量为【{cpu_count}】，大于该值，将被设置为{advise_cpu}。"
        info['并行'] = advise_cpu
    else:
        check['并行'] = 'OK'

    
    if info['线程'] > cpu_count:
        check['线程'] = f"CPU总核心数量为【{cpu_count}】，大于该值，将被设置为{advise_thread}。"
        info['线程'] = advise_thread
    else:
        check['线程'] = 'OK'

    if info['并行'] * info['线程'] > cpu_count:
        info['线程'] = advise_thread
        info['并行'] = advise_cpu
        check['并行'] = f"CPU总核心数量为【{cpu_count}】，【并行】的【线程】需要小于该值。"
        check['线程'] = f"CPU总核心数量为【{cpu_count}】，【并行】的【线程】需要小于该值。"
    
    return pd.DataFrame.from_dict(check,orient='index')[0]

def run_ana(info):
    # batch, outdir, baseline_list, parralel=2
    inpath = os.path.join(app.config['RAW_PATH'],info['批次'])
    outpath = os.path.join(app.config['OUT_PATH'], info['输出'])
    parralel = info['并行']
    thread =  info['线程']
    base = info['基线']
    baseline_list = []
    for b in base.split(';'):
        baseline_list.append('-bg')
        baseline_list.append(os.path.join(app.config['BASE_PATH'],str(b)+'.txt'))
    binpath = app.config['BIN_PATH']
    db = app.config['DB_PATH']
    # cmd = ['python', binpath, '-i', inpath, '-o', outpath]
    if os.path.exists(outpath):
        shutil.rmtree(outpath)
    os.makedirs(outpath)

    cmd = ['python', binpath, '-p', '_', '-i', inpath, '-o', outpath, 
    '-d', db, '-m', str(parralel), '-t', str(thread)]+baseline_list
    # '-bg', ' '.join(baseline_list)]
    print(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return p.pid

def log2html(logfile):
    log_dict = {}
    outdir = os.path.split(logfile)[0]
    with open(logfile,'rt') as h:
        for l in h.readlines():
            tmp = l.strip().split('\t')
            log_dict[tmp[0]] = log_dict.get(tmp[0],{})
            log_dict[tmp[0]][tmp[1]] = tmp[2]
    df = pd.DataFrame.from_dict(log_dict)#.loc[['run','DONE']]
    df_sta = pd.DataFrame()
    for s in os.listdir(outdir):
        stafile = os.path.join(os.path.join(outdir, s, 'sta.csv'))
        if os.path.isfile(stafile):
            df_tmp = pd.read_csv(stafile,header=None,index_col=0)
            df_tmp.columns = [s]
            df_sta = pd.concat([df_sta, df_tmp], axis=1)
            # if df_sta.empty:
                # df_sta = df_tmp
            # else:
                # df_sta = pd.merge(df_sta, df_tmp, left_index=True, right_index=True, how='outer')
    
    df = pd.concat([df, df_sta]).T.reset_index()
    # df = pd.read_csv(logfile,sep='\t',header=None,index_col=0)
    return df.to_html(na_rep='',)


cpu_count = int(multiprocessing.cpu_count())
advise_cpu = (cpu_count-2)/4
if advise_cpu < 1:
    advise_cpu = 1
else:
    advise_cpu = int(advise_cpu)
advise_thread = 4

@app.route('/getinfo', methods=['GET', 'POST'])
def getInfo():
    allBases = get_baseline()
    allBatch = get_batches()
    selectedChoices = ChoiceObj('bases', session.get('selected'))
    form = BasesLookupForm(obj=selectedChoices)
    form.bases.choices =  [(c, c) for c in allBases]
    default_str = '选择一个批次'
    form.batch.choices = [('0',f'{default_str:-^50}')] + [(c, c) for c in allBatch]
    form.batch.default = 0 
    if form.validate_on_submit():
        if form.submit.data:
            session['selected'] = form.bases.data
            info = {'批次':form.batch.data,'输出':form.outdir.data,'基线':';'.join(form.bases.data),'并行':form.parallel.data, '线程':form.thread.data}
            return redirect(url_for('ana', info=json.dumps(info)))

        if form.search.data:
            info = {'输出':form.outdir.data}
            df_info = pd.DataFrame.from_dict(info,orient='index')
            df_info.columns=['信息确认']
            return render_template('res.html',info=df_info.to_html(col_space=50), info_dict=info)

    form.parallel.data = advise_cpu
    form.thread.data = advise_thread
    return render_template('info.html',
                           form=form,
                           selected=session.get('selected',[]))


@app.route('/ana/<info>', methods=['GET', 'POST'])
def ana(info):
    df_info = pd.DataFrame.from_dict(json.loads(info),orient='index')
    df_info.columns=['信息确认']
    df_info['说明'] = check_info(df_info)
    ana_form = AnaForm()
    if ana_form.validate_on_submit():
        info_dict = df_info['信息确认'].to_dict()
        pid = run_ana(info_dict)
        print(f"PID for {info_dict} is {pid}")
        return render_template('res.html',info=df_info.to_html(col_space=50), info_dict=info_dict)
    return render_template('ana.html',info=df_info.to_html(col_space=50),form=ana_form)

@app.route('/get_ana_res', methods=['GET','POST'])
def get_ana_res():
    info = request.form.to_dict()
    outlog = os.path.join(app.config['OUT_PATH'], info['输出'], 'log')
    if os.path.isfile(outlog):
        df_html = log2html(outlog)
        # df = pd.read_csv(outlog,sep='\t')
        return df_html
    else:
        return f"【{info['输出']}】中【log】不存在"

@app.route('/test',methods=['GET','POST'])
def test():
    df_info = pd.DataFrame.from_dict({'a':'a','b':'b'},orient='index')
    form_ana = AnaForm()
    if form_ana.validate_on_submit():
        return ('test')

    return render_template('ana.html',info=df_info.to_html(col_space=20),form=form_ana)
