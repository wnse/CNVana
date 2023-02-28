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
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
# import numpy as np
plt.rcParams['font.sans-serif']=['SimHei']
# plt.rc("font",family='Arial')
# %%
def get_chr_num(c):
    n = c.lstrip('chr')
    if n in ['X','x']:
        n = 23
    if n in ['Y','y']:
        n = 24
    return int(n)


# %%
def plot_chr_scatter_rect(df, outpng):
    x_label = []
    x_label_pos = []
    x_label_color = []
    last_pos = 0
    color_list = ['#009900','#FFCC00','#990099','#3333CC']
    fig, ax = plt.subplots(figsize=(30,4))

    for i in df['chr_num'].sort_values().unique():
        df_tmp = df[df['chr_num']==i]
        df_tmp = df_tmp.sort_values(by='Position')
        df_tmp['Position'] = df_tmp['Position'] + last_pos
        last_pos = df_tmp['Position'].max()
        c = i
        if i == 1:
            c = 'Chr1'
        if i == 23:
            c = 'X'
        if i == 24:
            c = 'Y'
        x_label.append(c)
        x_label_pos.append(df_tmp['Position'].median())
        color = color_list[i%4-1]
        x_label_color.append(color)

        plt.scatter(df_tmp['Position'], df_tmp['copyNum'], s=15, c=color )

    ax.set_yticks(range(0,7,1))
    ax.set_yticklabels([0,'',2,'',4,'',6], fontdict={'size':20})
    ax.set_xticks(x_label_pos)
    ax.set_xticklabels(x_label, fontsize=20)
    for i, t in zip(x_label_color, ax.xaxis.get_ticklabels()):
        t.set_color(i)

    ax.set_ylabel('Copy Number',fontdict={'size':30})

    ax.set_xlim(0,last_pos+1e8)
    ax.tick_params(which='major',direction='in')
    ax.grid(visible=None, which='major', axis='y', color='grey', linestyle='dotted',linewidth=1)
    plt.savefig(outpng,bbox_inches='tight',pad_inches=0)
    

def plot_chr_scatter_rect_old(df, outpng):
    x_label = []
    x_label_pos = []
    x_label_color = []
    last_pos = 0
    color_list = ['#009900','#FFCC00','#990099','#3333CC']
    # fig, ax = plt.subplots(figsize=(30,4))
    fig = plt.figure(figsize=(30,4))
    ax = fig.add_axes([0.02,0.35,0.975,0.64])
    for i in df['chr_num'].sort_values().unique():
        df_tmp = df[df['chr_num']==i]
        df_tmp = df_tmp.sort_values(by='Position')
        df_tmp['Position'] = df_tmp['Position'] + last_pos
        last_pos = df_tmp['Position'].max()
        c = i
        if i == 1:
            c = 'Chr1'
        if i == 23:
            c = 'X'
        if i == 24:
            c = 'Y'
        x_label.append(c)
        x_label_pos.append(df_tmp['Position'].median())
        color = color_list[i%4-1]
        x_label_color.append(color)

        plt.scatter(df_tmp['Position'], df_tmp['copyNum'], s=5, c=color )

    X_l=[0.1,0.3,0.5,0.6,0.75,0.9]
    Y_l=[0.65,0.55,0.8,0.5,0.85,0.45]
    for x_tmp,y_tmp in zip(X_l, Y_l):
        # ax.annotate(text="仅供科学研究", xy=(1e8, 4.5), fontsize=40, color='grey',style='italic')
        ax.text(x=x_tmp-0.07, y=y_tmp-0.3, s="仅供科学研究", fontsize=40, color='grey',transform=ax.transAxes,
            fontstyle='italic',fontweight="bold",alpha=0.3)

    ax.set_yticks(range(0,7,1))
    ax.set_yticklabels([0,'',2,'',4,'',6], fontdict={'size':20})
    ax.set_xticks(x_label_pos)
    ax.set_xticklabels(x_label, fontsize=20)
    for i, t in zip(x_label_color, ax.xaxis.get_ticklabels()):
        t.set_color(i)

    ax.set_ylabel('Copy Number',fontdict={'size':30})

    ax.set_xlim(0,last_pos+5e7)
    ax.set_ylim(0,6.2)
    ax.tick_params(which='major',direction='in')
    ax.grid(visible=None, which='major', axis='y', color='grey', linestyle='dotted',linewidth=1)
    plt.savefig(outpng,dpi=300)#,bbox_inches='tight',pad_inches=0)
    

# %%
def plot_chr_scatter_square(df, outpng, color_num=4):
    color_list = ['#009900','#FFCC00','#990099','#3333CC']        
    if color_num==2:
        color_list = ['red','blue']
    fig, ax = plt.subplots(3,1,figsize=(40,18))
    
    def plot_chr(df_chr, ax):
        x_label = []
        x_label_pos = []
        x_label_color = []
        last_pos = 0
        for i in df_chr['chr_num'].sort_values().unique():
            df_tmp = df[df['chr_num']==i]
            df_tmp = df_tmp.sort_values(by='Position')
            x_label.extend(range(0,df_tmp['Position'].max(),30000000)[1:])
            df_tmp['Position'] = df_tmp['Position'] + last_pos
            x_label_pos.extend(range(last_pos,df_tmp['Position'].max(),30000000)[1:])
            color = color_list[i%len(color_list)-1]        

            c = i
            if i == 23:
                c = 'X'
            if i == 24:
                c = 'Y'
            c = 'chr'+str(c)
            ax.scatter(df_tmp['Position'], df_tmp['copyNum'], s=15, c=color )
            ax.annotate(text=c, xy=(df_tmp['Position'].min(), 4.5), fontsize=20, color=color)
            ax.plot([last_pos, df_tmp['Position'].max()],[0,0],'-',color=color, linewidth=3)
            last_pos = df_tmp['Position'].max()

        ax.set_yticks(range(0,5,1))
        ax.set_yticklabels([str(i)+'X' for i in range(0,5,1)], fontdict={'size':15})
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_position(('data',0))
        ax.plot([0,0],[0,5],'-',color='red',linewidth=3)
    #     ax.tick_params(axis='y', colors='red')
        ax.set_xticks(x_label_pos)
        ax.set_xticklabels([str(int(i/1e6)) for i in x_label], fontsize=15)
        for i, t in zip(x_label_color, ax.xaxis.get_ticklabels()):
            t.set_color(i)

    #     ax.set_ylabel('Copy Number',fontdict={'size':40})

        ax.set_xlim(0,last_pos+1e7)
        ax.tick_params(which='major',direction='in')
        ax.grid(visible=None, which='major', axis='y', color='grey', linestyle='dotted',linewidth=1)
#     return ax

    plot_chr(df[df['chr_num'].isin(range(1,6))], ax[0])
    plot_chr(df[df['chr_num'].isin(range(6,13))], ax[1])
    plot_chr(df[df['chr_num'].isin(range(13,25))], ax[2])

    plt.subplots_adjust(hspace=0.01)
#     plt.show()
    plt.savefig(outpng,bbox_inches='tight',pad_inches=0) 


# %%
def get_plot(csv_file, plot_type=1, png_file=None, sex_chr=''):
    if not png_file:
        png_file = csv_file + '.png'
#     csv_file  = '/Users/yangkai/work/GermountX/CNV/广州/2023_01_17/B_PGTA_ZY06.fq_merge_bin_correct_da.csv'
    df = pd.read_csv(csv_file)
    df['chr_num'] = df['chr'].apply(get_chr_num)
    if sex_chr == 'XX':
        df.loc[df['chr'].isin(['y','Y','chrY']), 'copyNum']=0
    df = df.sort_values(by=['chr_num', 'Position']).reset_index()
    if plot_type==0:
        plot_chr_scatter_rect_old(df, png_file)
    if plot_type==1:
        plot_chr_scatter_rect(df, png_file)
    if plot_type==2:
        plot_chr_scatter_square(df, png_file)
    if plot_type==3:
        plot_chr_scatter_square(df, png_file, color_num=2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='cnv file (.csv)')
    parser.add_argument('-t', '--type', default=1, type=int, choices=[0,1,2,3], help='plot type for rect or squre or color')
    parser.add_argument('-o', '--output', default=None, help='out png file, default input.png')
    parser.add_argument('-s', '--sex', default='', type=str, help='plot sex chr')
    args = parser.parse_args()
    get_plot(args.input, args.type, args.output, args.sex)

# %%
# csv_file  = '/Users/yangkai/work/GermountX/CNV/广州/2023_01_17/B_PGTA_ZY06.fq_merge_bin_correct_da.csv'
# df = pd.read_csv(csv_file)

# df['chr_num'] = df['chr'].apply(get_chr_num)
# df = df.sort_values(by=['chr_num', 'Position']).reset_index()

# # plot_chr_scatter_width(df, 'test_out.png')
# plot_chr_scatter_square(df, 'test_square.png',color_num=2)

