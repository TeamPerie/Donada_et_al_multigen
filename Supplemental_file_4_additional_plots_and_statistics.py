##########################################################################################
### Overview: This script is to generate some exploratory plots from the heatmap, and to perform               statistical comparisons between experimental conditions
###
### Correspondence: Leila Perie (leila.perie@curie.fr) Alessandro Donada (alessandro.donada@gmail.com)
###
### Date: 09/12/2022
##################################################################


#!/usr/bin/env python
# coding: utf-8

# In[1]: 2_bar_plot

# Step  one:  prepare the  workspace, defining all libraries, plotting parameters, and paths that we  will use in our analysis

#import all of the libraries that are needed for this  analysis

import copy as cp
import itertools as itt
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
import pandas as pd
import pickle
import scipy.optimize as spo
import scipy.special as spsp
import scipy.stats as sps
import seaborn as sns
from matplotlib.ticker import FixedLocator

# set the plotting parameters for the matplotlib library
mpl.rcParams['axes.titlesize'] = 'xx-large'
mpl.rcParams['axes.labelsize'] = 'xx-large'
mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.framealpha'] = 0.5
mpl.rcParams['legend.fontsize'] = 'large'


#set the pathways for input (dataset tables) and output (modified tables and figures)
mydir = os.path.dirname(os.getcwd())
path_proj = os.path.join(mydir, 'MultiGen_analysis')
path_sc = os.path.join(path_proj, 'csv', 'Single_cell')
path_pickle = os.path.join(path_proj, 'pickled_data')
path_sort = os.path.join(path_proj, 'csv', 'Sort')
path_plot = os.path.join(path_proj, 'figures')
print(path_proj)



# In[2]:

#Load the processed data

df = pd.read_csv(os.path.join(path_sc, 'Pooled_data.csv'), sep=',', decimal='.')


# In[3]:

#This cell set the colours for the different cell types and generations, organising the specific legends. cell_cols refer to the
#cell types detected and specified in the heatmap. gen_cols refer to the number of divisions (generations) measured experimentally,
respectively from generation 0 (no division) to generation 6. cnd_cols is referring to eventual columns describing two experimental
conditions. div1_cols is specifically referring to the type of division performed during division 1, specified in sym_labs. 

#COLORS
cell_cols = {'HSC':'seagreen', 'MPP':'Limegreen', 'LMPP':'Red', 'CMP':'Deepskyblue','GMP':'Violet', 'MEP':'Royalblue', 'CD34-':'aqua'}
gen_cols = np.array(["#660000", "#FF1000", "#FF9F00", "#40FFBF", "#00AFFF", "#0020FF",
                     "#000078"])
cnd_cols = np.array(['Steelblue', 'Gold'])
div1_cols = np.array(['Royalblue', 'paleturquoise', 'Gold', 'khaki'])


# In[4]:

#(EXPERIMENT SPECIFIC)        
#This cell set few experiment specific parameters and organise the display of the cell types

#Cell class supports and colours (dictionary values) ordered by progenitors (dictionary keys) 
class_dct = {'HSC':['HSC','MPP','LMPP','CMP','GMP','MEP','CD34-'],
           'MPP':['HSC','MPP','LMPP','CMP','GMP','MEP','CD34-'],
           'HPC':['HSC','MPP','LMPP','CMP','GMP','MEP','CD34-']}

temp_col = np.array(list(cell_cols.values()))
#temp_col = np.array(temp_col[0:1]+temp_col[2:])
cols_dct = {'HSC':temp_col,
            'MPP':temp_col[[0,1,2,3,4,5,6]],
            'HPC':temp_col[[0,1,2,3,4,5,6]]}
            
#Iterating/data filtering lists. Those filters must be specified by the user: conds stands for the experimental conditions, or_cells
#to the cell of origin, times to the cell culture length. sym_labs refers to the type of division 1, specified further down.

conds = ['GT', 'Diff'] 
or_cells = np.unique(df.Original_cell)
sym_labs = ['Sym Undiff', 'Sym Diff', 'Asym Undiff', 'Asym Diff']
times = ['72h']
or_cells


# In[5]:

#Function that calculate a statistical function (func), with fixed arguments (args) on the actual data.
#For each element of the statistic, the confidence intervals are calculated by basic bootstrapping, using a number of iterations (n_iter)
#IF bool_trim01=TRUE, the confidence intervals are capped within then [0,1] interval

def stat_bootCI(func, args, data, n_iter=250000, bool_trim01=True):
    stat = func(data, *args)
    boot_stat = np.array([func(data[boot_index], *args)
                          for boot_index in np.random.choice(len(data), size=(n_iter,len(data)))])
    boot_per = np.percentile(boot_stat, q=[97.5, 2.5], axis=0)
    boot_ci = 2.*stat - boot_per
    if bool_trim01:
        boot_ci[0,boot_ci[0,:]<0.] = np.zeros((boot_ci[0,:]<0.).sum())
        boot_ci[1,boot_ci[1,:]>1.] = np.ones((boot_ci[1,:]>1.).sum())
    return stat, boot_ci

#Statistics for FIG 2B-2E
#CUMULATIVE PROPORTION OF CELL CLASSES ORDERLY IN ordered_classes
def class_distribution(data, ordered_classes, bool_stack=True):
    if bool_stack:
        classes, counts = np.unique(np.hstack(data), return_counts=True)
    else:
        classes, counts = np.unique(data, return_counts=True)
    counts_sum = counts.sum()
    if counts_sum == 0:
        return np.zeros(len(ordered_classes))
    else:
        return np.cumsum([counts[classes==cl] if any(classes==cl) else np.zeros(1) for cl in ordered_classes]
                    )/counts_sum



### STEP 1: produce the parameters (dct_plot) for the bar plot in each panel, then pickle it for STEP 2.


# In[6]:

#FIG 2B: aggregate frequency per cell type.
        
dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (class_dct[oc], True)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)]
        data = np.array([df_temp[df_temp.Family==fam].Class.values for fam in np.unique(df_temp.Family)])
        print(data)
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'Fig2B'), 'wb') as fp:
    pickle.dump(dct_plot, fp)


# In[7]:

#FIG 2C: aggregate percentage of families per generation

dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (np.arange(7), False)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)]
        data = np.array([df_temp[df_temp.Family==fam].Generation.values.max() for fam in np.unique(df_temp.Family)])
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'Fig2C'), 'wb') as fp:
    pickle.dump(dct_plot, fp)


# In[8]:

#FIG 2D: frequency of cell type per generation 1, same as Fig 2B but restricted to cells in generation 1
    
dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (class_dct[oc], True)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)&(df.Generation==1)]
        data = np.array([df_temp[df_temp.Family==fam].Class.values for fam in np.unique(df_temp.Family)])
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'SuppFig2'), 'wb') as fp:
    pickle.dump(dct_plot, fp)
    
    
# In[9]:

#FIG 2E: label a family (df DataFrame), assumed having 2 cells (daughter 1 and daughter 2) both in Generation 1. 
#The divisions patterns are specified with sym_labs, and refer to the type of division performed. The first_div_class function
#use a loop to verify the type of division performed: if both cells have the same phenotype and it is the same phenotype of the 
#cell of origin (or_cell), means that the division is a symmetric undifferentiated division. If both cells have the same phenotype, 
#and it is not the same phenotype of the cell of origin (or_cell), means that the division is a symmetric differentiated division.
#If the two daughter cells have a different phenotype, but at least one of them share the same phenotype of the mother, the division
#is an asymmetric undifferentiated division. Finally, if both daughter cells have different phenotypes between themselves and the 
#mother cells, means that the division is asymmetric differentiated.

def first_div_class(df):
    daugher1, daugher2 = df.Class.values
    progenitor = df.Original_cell.values[0]
    if daugher1 == daugher2:
        if daugher1 == progenitor:
            return sym_labs[0]
        else:
            return sym_labs[1]
    else:
        if (daugher1 == progenitor) or (daugher2 == progenitor):
            return sym_labs[2]
        else:
            return sym_labs[3]

dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (sym_labs, False)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)&(df.Generation==1)]
        data = np.array([first_div_class(df_temp[df_temp.Family==fam])
                         for fam in np.unique(df_temp.Family) if len(df_temp[df_temp.Family==fam])==2])
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'Fig2E'), 'wb') as fp:
    pickle.dump(dct_plot, fp)


# In[10]:

#FIG 3A: 
#STATISTIC FOR FIG 3A, SAME AS class_distribution,
#BUT CALCULATED ON CELLS IN GENERATIONS 0 TO 3 (generations[:-1]) AND 4+ (generations[-1]) separately

def class_distribution_gen(data, ordered_classes, generations):
    data_stack = np.vstack(data)
    res = [class_distribution(data_stack[data_stack[:,1]==g,0], ordered_classes, bool_stack=False)
           for g in generations[:-1]]
    res = np.array(
        res + [class_distribution(data_stack[data_stack[:,1]>=generations[-1],0],
                                  ordered_classes, bool_stack=False)]
    )
    return res

generations = np.arange(5, dtype=int)


lst_oc_t_cnd = [(oc,t,cnd) for oc in np.unique(df.Original_cell)
                   for t in np.unique(df.Culture_time)
                   for cnd in np.unique(df.Culture_condition)]

for oc,t,cnd in lst_oc_t_cnd:
    dct_plot = {}
    args = (class_dct[oc], generations)
    df_temp = df[(df.Original_cell==oc)&(df.Culture_time==t)&(df.Culture_condition==cnd)]
    data = np.array([df_temp[df_temp.Family==fam].loc[:,['Class','Generation']].values
                     for fam in np.unique(df_temp.Family)])
    stat, ci = stat_bootCI(class_distribution_gen, args, data)    
    for gen in generations:
        dct_plot[str(gen)] = {'':(stat[gen], ci[:,gen],)}
    with open(os.path.join(path_pickle, 'Fig3A' +oc+t+cnd), 'wb') as fp:
        pickle.dump(dct_plot, fp)



# In[11]:

#FIG 3B: frequency of cell type per generation 0, same as Fig 2B but restricted to cells in generation 0

dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (class_dct[oc], False)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)&(df.Generation==0)]
        data = np.array([df_temp[df_temp.Family==fam].Class.values for fam in np.unique(df_temp.Family)])
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'Fig3B'), 'wb') as fp:
    pickle.dump(dct_plot, fp)

# In[12]:

#FIG 3C: frequency of cell type per generation 2, same as Fig 2B but restricted to cells in generation 2

dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (class_dct[oc], True)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)&(df.Generation==2)]
        data = np.array([df_temp[df_temp.Family==fam].Class.values for fam in np.unique(df_temp.Family)])
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'Fig3C'), 'wb') as fp:
    pickle.dump(dct_plot, fp)



# In[13]:

#FIG 3D: frequency of cell type per generation 3, same as Fig 2B but restricted to cells in generation 3

dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (class_dct[oc], True)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)&(df.Generation==3)]
        data = np.array([df_temp[df_temp.Family==fam].Class.values for fam in np.unique(df_temp.Family)])
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'Fig3D'), 'wb') as fp:
    pickle.dump(dct_plot, fp)



# In[14]:

#FIG 3E: frequency of cell type per generation 4, same as Fig 2B but restricted to cells in generation 4

dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (class_dct[oc], True)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)&(df.Generation==4)]
        data = np.array([df_temp[df_temp.Family==fam].Class.values for fam in np.unique(df_temp.Family)])
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'Fig3E'), 'wb') as fp:
    pickle.dump(dct_plot, fp)


# In[15]:

#FIG 3F: frequency of cell type per generation 5, same as Fig 2B but restricted to cells in generation 5

dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (class_dct[oc], True)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)&(df.Generation==5)]
        data = np.array([df_temp[df_temp.Family==fam].Class.values for fam in np.unique(df_temp.Family)])
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'Fig3F'), 'wb') as fp:
    pickle.dump(dct_plot, fp)


# In[16]:

#FIG 3G: frequency of cell type per generation 6, same as Fig 2B but restricted to cells in generation 6

dct_plot = {}
for oc in or_cells:
    dct_plot[oc] = {}
    args = (class_dct[oc], True)
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)&(df.Generation==6)]
        data = np.array([df_temp[df_temp.Family==fam].Class.values for fam in np.unique(df_temp.Family)])
        dct_plot[oc][cnd] = stat_bootCI(class_distribution, args, data)
with open(os.path.join(path_pickle, 'Fig3G'), 'wb') as fp:
    pickle.dump(dct_plot, fp)
    

# In[17]:
    
### STEP 2: plot the pickled data individually. Each block will code for the plot and the corresponding statistical test

def plot_bars(dct_plot, cols_dct, xlabel, ylabel,
              legend_title, legend_classes, legend_cols,
              cat_width=0.7, str_sep='   '):
    fig, ax = plt.subplots()
    n_cat = len(dct_plot.keys())
    for nc,cat_k in enumerate(dct_plot.keys()):
        x_cat = -cat_width/2 + nc
        n_hue = len(dct_plot[cat_k].keys())
        for nh, hue_k in enumerate(dct_plot[cat_k].keys()):
            y = dct_plot[cat_k][hue_k][0]
            ci = dct_plot[cat_k][hue_k][1]
            bar_width = cat_width/n_hue
            bottoms = np.hstack([np.zeros(1),y[:-1]])
            x_hue = x_cat + nh*bar_width
            x_ci = [x_hue + k*bar_width/(len(y)) for k in np.arange(1,len(y))]
            if len(ci) != 0:
                for k in np.arange(len(y)-1):
                    ax.plot([x_ci[k]]*3, [ci[0,k], y[k], ci[1,k]],
                            ls='-', marker='_', ms=5., lw=1., color='Black')
            ax.bar(x=[x_hue]*len(y), height=y-bottoms,
                   width=[bar_width]*len(y), bottom=bottoms, align='edge',
                   color=cols_dct[cat_k])
    handles = [mpatches.Patch(color=co, label=cl) for cl,co in zip(legend_classes,legend_cols)]
    ax.legend(handles=handles, loc='upper left', ncol=1, bbox_to_anchor= (1.05, 1.), title=legend_title)

    ax.set_xticks(np.arange(n_cat, dtype=int))
    xtickslabels = [str_sep.join(dct_plot[cat_k].keys())+'\n'+cat_k for cat_k in dct_plot.keys()]
    ax.set_xticklabels(xtickslabels)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_yticks(np.linspace(0,1,5))
    ax.set_yticklabels(np.linspace(0,100,5, dtype=int))
    sns.despine(fig=fig)
    return fig,ax


# In[18]:

#FIG 2B: aggregate frequency per cell type
            
with open(os.path.join(path_pickle,'Fig2B'), 'rb') as fp:
    dct_plot = pickle.load(fp)
with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, cols_dct, xlabel='Progenitor type', ylabel='Cells (%)', legend_title='Cell type',
                       legend_classes=class_dct['HSC'], legend_cols=cols_dct['HSC'])
    ax.set_title('Cell type percentages')
    fig.savefig(os.path.join(path_plot,'Cell type percentages (total).pdf'), bbox_inches='tight')

print('New TEST DIFFERENTIATION DEPENDENCE ON CONDITION UNCORR, now using G-test stat')
def spread_support(support, data):
    res = np.zeros(len(support))
    d_sup,freq = np.unique(data, return_counts=True)
    res[[k in d_sup for k in support]] = freq
    return res
mc_iter = 250000.
#dft, dfi = dct_data['pooled']
df_temp = cp.deepcopy(df)

for n_oc,oc in enumerate(or_cells):
    df_cnd = df_temp[(df_temp.Original_cell==oc)]
    df_ary = np.array(df_cnd.loc[:,['Class','Family','Culture_condition']])
    cls_lst = np.unique(df_ary[:,0])
    fams_lst = np.unique(df_ary[:,1])
    fam_cl = np.array([df_ary[df_ary[:,1]==fam,0] for fam in fams_lst])
    fam_cnd = np.array([df_ary[df_ary[:,1]==fam,2][0] for fam in fams_lst])
    fam_data = np.array(list(zip(fam_cl, fam_cnd)))
    support = np.unique(df_ary[:,0])
    table = np.array([
            spread_support(support, np.hstack(fam_data[fam_data[:,1]==cnd,0]))
            for cnd in conds])
    stat = sps.chi2_contingency(table, correction=False, lambda_='log-likelihood')[0]
    per_indx = [np.random.permutation(fam_data[:,1]) for k in np.arange(mc_iter)]
    per_tables = np.array([
                [spread_support(support, np.hstack(fam_data[p_ind==cnd,0])) for cnd in conds]
            for p_ind in per_indx])
    per_stat = np.array([sps.chi2_contingency(p_tab, correction=False, lambda_='log-likelihood')[0]
            for p_tab in per_tables])
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_'+oc, 'rb') as fp:
    #with open('.pickled_data/Fig2B', 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    fig,ax=plt.subplots()
    sns.distplot(per_stat, norm_hist=True, kde=False, ax=ax)
    ax.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    ax.set_xlabel('G-test statistic')
    ax.set_ylabel('Distribution')
    ax.text(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(fams_lst)),
                horizontalalignment='right', verticalalignment='top', fontsize='large',
                transform=fig.transFigure)
    ax.set_title('Test cell type percentages '+oc+' (total)')
    sns.despine(fig=fig)
    fig.savefig(os.path.join(path_plot,'Test cell type percentages '+oc+' (total).pdf'))


# In[19]:

#FIG 2C: aggregate percentage of families per generation

with open(os.path.join(path_pickle,'Fig2C'), 'rb') as fp:
    dct_plot = pickle.load(fp)
gen_cols_dct = {k:gen_cols for k in dct_plot.keys()}

with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, gen_cols_dct, xlabel='Progenitor', ylabel='Families (%)', legend_title='Generation',
                       legend_classes=[str(gen) for gen in np.arange(7, dtype=int)], legend_cols=gen_cols)
    ax.set_title('Family percentages per generations')
    fig.savefig(os.path.join(path_plot,'Family percentages per generations.pdf'), bbox_inches='tight')
    
print('TEST GENERATION PROGRESSION DEPENDENCE ON CONDITION')
mc_iter = 250000.
def cnd_stat(lab_cnd, sym_labs, conds, k=0.):
    table = np.array([[(lab_cnd[lab_cnd[:,1]==cnd,0] == lab).sum() for lab in sym_labs]
              for cnd in conds])
    return sps.chi2_contingency(table, correction=False, lambda_="log-likelihood")[0]
for n_oc,oc in enumerate(or_cells):
    df_cnd = df[(df.Original_cell==oc)]
    lab_cnd = np.array([(df_cnd[df_cnd.Family==fam].Generation.max(), df_cnd[df_cnd.Family==fam].Culture_condition.iloc[0])
                        for fam in np.unique(df_cnd.Family)])
    sym_labs = [str(g) for g in np.arange(df_cnd.Generation.max()+1) if str(g) in lab_cnd[:,0]]
    print(sym_labs)
    stat = cnd_stat(lab_cnd, sym_labs, conds)
    per_stat = np.array([cnd_stat(np.vstack((lab_cnd[:,0], np.random.permutation(lab_cnd[:,1]))).T,
                          sym_labs, conds, k) for k in np.arange(mc_iter)])
    with open('./pickled_data/Cnd_maxdd_condition_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/Cnd_maxdd_condition_'+oc, 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    plt.figure()
    sns.distplot(per_stat, norm_hist=True, kde=False)
    plt.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    plt.xlabel('Log-likelihood ratio statistic')
    plt.ylabel('Distribution')
    plt.figtext(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(lab_cnd)),
                horizontalalignment='right', verticalalignment='top', fontsize='large')
    plt.title('Family percentages per generations '+oc)
    sns.despine()
    plt.savefig(os.path.join(path_plot,'Test family percentages per generations '+oc+'.pdf'))
    plt.show()
    plt.close()    


# In[20]:
    
#FIG 2D: frequency of cell type per generation 1

with open(os.path.join(path_pickle,'SuppFig2'), 'rb') as fp:
    dct_plot = pickle.load(fp)
with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, cols_dct, xlabel='Progenitor', ylabel='Cells (%)', legend_title='Cell type',
                       legend_classes=class_dct['HSC'], legend_cols=cols_dct['HSC'])
    ax.set_title('Cell types percentages Generation 1')
    fig.savefig(os.path.join(path_plot,'Cell type percentages (Gen1).pdf'), bbox_inches='tight')
    

def spread_support(support, data):
    res = np.zeros(len(support))
    d_sup,freq = np.unique(data, return_counts=True)
    res[[k in d_sup for k in support]] = freq
    return res
mc_iter = 250000.
for n_oc,oc in enumerate(or_cells):
    df_cnd = df_temp[(df_temp.Original_cell==oc)&(df_temp.Generation==1)]
    df_ary = np.array(df_cnd.loc[:,['Class','Family','Culture_condition']])
    cls_lst = np.unique(df_ary[:,0])
    fams_lst = np.unique(df_ary[:,1])
    fam_cl = np.array([df_ary[df_ary[:,1]==fam,0] for fam in fams_lst])
    fam_cnd = np.array([df_ary[df_ary[:,1]==fam,2][0] for fam in fams_lst])
    fam_data = np.array(list(zip(fam_cl, fam_cnd)))
    support = np.unique(df_ary[:,0])
    table = np.array([
            spread_support(support, np.hstack(fam_data[fam_data[:,1]==cnd,0]))
            for cnd in conds])
    stat = sps.chi2_contingency(table, correction=False, lambda_='log-likelihood')[0]
    per_indx = [np.random.permutation(fam_data[:,1]) for k in np.arange(mc_iter)]
    per_tables = np.array([
                [spread_support(support, np.hstack(fam_data[p_ind==cnd,0])) for cnd in conds]
            for p_ind in per_indx])
    per_stat = np.array([sps.chi2_contingency(p_tab, correction=False, lambda_='log-likelihood')[0]
            for p_tab in per_tables])
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen1_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen1_'+oc, 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    fig,ax=plt.subplots()
    sns.distplot(per_stat, norm_hist=True, kde=False, ax=ax)
    ax.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    ax.set_xlabel('G-test statistic')
    ax.set_ylabel('Distribution')
    ax.text(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(fams_lst)),
                horizontalalignment='right', verticalalignment='top', fontsize='large',
                transform=fig.transFigure)
    ax.set_title('Test cell type percentages Generation 1 '+oc)
    sns.despine(fig=fig)
    fig.savefig(os.path.join(path_plot,'Test cell type percentages '+oc+' (Gen1).pdf'))


# In[21]:
    
#FIG 2E: LABEL A FAMILY (df DataFrame, ASSUMED HAVING 2 CELLS daugher1 AND daugher2 IN GENERATION 1)
#AFTER ONE OF THE DIVISION PATTERNS IN sym_labs    

with open(os.path.join(path_pickle,'Fig2E'), 'rb') as fp:
    dct_plot = pickle.load(fp)
div1_cols_dct = {k:div1_cols for k in dct_plot.keys()}

with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, div1_cols_dct, xlabel='Progenitor', ylabel='Generation 1 families (%)',
                       legend_title='First division type', legend_classes=sym_labs,
                       legend_cols=div1_cols_dct['HSC'])
    ax.set_title('Division type Generation 1')
    fig.savefig(os.path.join(path_plot,'Division type (Gen1).pdf'), bbox_inches='tight')
   
    
print('TEST GEN 1 (A)SYMMETRY DEPENDENCE ON CONDITION')
mc_iter = 250000.
def cnd_stat(lab_cnd, sym_labs, conds, k=0.):
    table = np.array([[(lab_cnd[lab_cnd[:,1]==cnd,0] == lab).sum() for lab in sym_labs]
              for cnd in conds])
    return sps.chi2_contingency(table, correction=False, lambda_="log-likelihood")[0]

def check_sym_asym(class_lst, progenitor):
    daugher1 = class_lst[0]
    daugher2 = class_lst[1]
    if daugher1 == daugher2:
        if daugher1 == progenitor:
            return sym_labs[0]
        else:
            return sym_labs[1]
    else:
        if (daugher1 == progenitor) or (daugher2 == progenitor):
            return sym_labs[2]
        else:
            return sym_labs[3]

for n_oc,oc in enumerate(or_cells):
    df_cnd = df_temp[(df_temp.Original_cell==oc)]
    lab_cnd = np.array([(check_sym_asym(df_cnd[df_cnd.Family==fam].Class.tolist(), oc),
                         df_cnd[df_cnd.Family==fam].Culture_condition.iloc[0])
                        for fam in np.unique(df_cnd.Family)
                            if len(df_cnd[df_cnd.Family==fam])==2])
    sym_labs_restr = [lab for lab in sym_labs if (lab in lab_cnd[:,0])]
    stat = cnd_stat(lab_cnd, sym_labs_restr, conds)
    per_stat = np.array([cnd_stat(np.vstack((lab_cnd[:,0], np.random.permutation(lab_cnd[:,1]))).T,
                          sym_labs_restr, conds, k) for k in np.arange(mc_iter)])
    with open('./pickled_data/gen1_condition_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/gen1_condition_'+oc, 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    
    plt.figure()
    sns.distplot(per_stat, norm_hist=True, kde=False)
    plt.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    plt.xlabel('Log-likelihood ratio statistic')
    plt.ylabel('Distribution')
    plt.figtext(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(lab_cnd)),
                horizontalalignment='right', verticalalignment='top', fontsize='large')
    plt.title('Test Division type Gen 1 '+oc)
    sns.despine()
    plt.savefig(os.path.join(path_plot,'Test Division type '+oc+' Gen1.pdf'))
    plt.show()
    plt.close()    


# In[22]:
    
#FIG 3A: frequency of cell type per all the generations detected in the experiment. 
#NOTE: as not all generations would be present, check beforehand which generations are missing.
    
for oc,t,cnd in lst_oc_t_cnd:
    with open(os.path.join(path_pickle,'Fig3A'+oc+t+cnd), 'rb') as fp:
        dct_plot = pickle.load(fp)
    cols_gen_dct = {str(gen):cols_dct[oc] for gen in np.arange(7, dtype=int)}
    with sns.axes_style("ticks"):
        fig,ax = plot_bars(dct_plot, cols_gen_dct, xlabel='Generation', ylabel='Cells (%)',
                           legend_title='Cell type', legend_classes=class_dct['HSC'],
                           legend_cols=cols_dct['HSC'], cat_width=1., str_sep='')
        ax.set_title("Cell type percentages per generation "+cnd)
        fig.savefig(os.path.join(path_plot,'Cell type percentages per generation '+oc+' '+t+' '+cnd+'.pdf'), bbox_inches='tight')


# In[23]:
    
#FIG 3B: frequency of cell type per generation 0

with open(os.path.join(path_pickle,'Fig3B'), 'rb') as fp:
    dct_plot = pickle.load(fp)
with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, cols_dct, xlabel='Progenitor', ylabel='Cells (%)',
                       legend_title='Cell type', legend_classes=class_dct['HSC'],
                       legend_cols=cols_dct['HSC'])
    ax.set_title('Cell types percentages Generation 0')
    fig.savefig(os.path.join(path_plot,'Cell type percentages (Gen0).pdf'), bbox_inches='tight')

print('TEST GEN 0 CLASSES DEPENDENCE ON CONDITION')
def spread_support(support, data):
    res = np.zeros(len(support))
    d_sup,freq = np.unique(data, return_counts=True)
    res[[k in d_sup for k in support]] = freq
    return res

mc_iter = 250000.
for n_oc,oc in enumerate(or_cells):
    df_cnd = df_temp[(df_temp.Original_cell==oc)&(df_temp.Generation==0)]
    df_ary = np.array(df_cnd.loc[:,['Class','Family','Culture_condition']])
    cls_lst = np.unique(df_ary[:,0])
    fams_lst = np.unique(df_ary[:,1])
    fam_cl = np.array([df_ary[df_ary[:,1]==fam,0] for fam in fams_lst])
    fam_cnd = np.array([df_ary[df_ary[:,1]==fam,2][0] for fam in fams_lst])
    fam_data = np.array(list(zip(fam_cl, fam_cnd)))
    support = np.unique(df_ary[:,0])
    table = np.array([
            spread_support(support, np.hstack(fam_data[fam_data[:,1]==cnd,0]))
            for cnd in conds])
    stat = sps.chi2_contingency(table, correction=False, lambda_='log-likelihood')[0]
    per_indx = [np.random.permutation(fam_data[:,1]) for k in np.arange(mc_iter)]
    per_tables = np.array([
                [spread_support(support, np.hstack(fam_data[p_ind==cnd,0])) for cnd in conds]
            for p_ind in per_indx])
    per_stat = np.array([sps.chi2_contingency(p_tab, correction=False, lambda_='log-likelihood')[0]
            for p_tab in per_tables])
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen1_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen1_'+oc, 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    fig,ax=plt.subplots()
    sns.distplot(per_stat, norm_hist=True, kde=False, ax=ax)
    ax.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    ax.set_xlabel('G-test statistic')
    ax.set_ylabel('Distribution')
    ax.text(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(fams_lst)),
                horizontalalignment='right', verticalalignment='top', fontsize='large',
                transform=fig.transFigure)
    ax.set_title('Test cell type percentages '+oc+" gen0")
    sns.despine(fig=fig)
    fig.savefig(os.path.join(path_plot,'Test cell type percentages '+oc+' gen0.pdf'))
    

# In[24]:
    
#FIG 3C: frequency of cell type per generation 2

with open(os.path.join(path_pickle,'Fig3C'), 'rb') as fp:
    dct_plot = pickle.load(fp)
with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, cols_dct, xlabel='Progenitor', ylabel='Cells (%)',
                       legend_title='Cell type', legend_classes=class_dct['HSC'],
                       legend_cols=cols_dct['HSC'])
    ax.set_title('Cell types percentages Generation 2')
    fig.savefig(os.path.join(path_plot,'Cell type percentages (Gen2).pdf'), bbox_inches='tight')

print('TEST GEN 2 CLASSES DEPENDENCE ON CONDITION')
def spread_support(support, data):
    res = np.zeros(len(support))
    d_sup,freq = np.unique(data, return_counts=True)
    res[[k in d_sup for k in support]] = freq
    return res

mc_iter = 250000.
for n_oc,oc in enumerate(or_cells):
    df_cnd = df_temp[(df_temp.Original_cell==oc)&(df_temp.Generation==2)]
    df_ary = np.array(df_cnd.loc[:,['Class','Family','Culture_condition']])
    cls_lst = np.unique(df_ary[:,0])
    fams_lst = np.unique(df_ary[:,1])
    fam_cl = np.array([df_ary[df_ary[:,1]==fam,0] for fam in fams_lst])
    fam_cnd = np.array([df_ary[df_ary[:,1]==fam,2][0] for fam in fams_lst])
    fam_data = np.array(list(zip(fam_cl, fam_cnd)))
    support = np.unique(df_ary[:,0])
    table = np.array([
            spread_support(support, np.hstack(fam_data[fam_data[:,1]==cnd,0]))
            for cnd in conds])
    stat = sps.chi2_contingency(table, correction=False, lambda_='log-likelihood')[0]
    per_indx = [np.random.permutation(fam_data[:,1]) for k in np.arange(mc_iter)]
    per_tables = np.array([
                [spread_support(support, np.hstack(fam_data[p_ind==cnd,0])) for cnd in conds]
            for p_ind in per_indx])
    per_stat = np.array([sps.chi2_contingency(p_tab, correction=False, lambda_='log-likelihood')[0]
            for p_tab in per_tables])
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    fig,ax=plt.subplots()
    sns.distplot(per_stat, norm_hist=True, kde=False, ax=ax)
    ax.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    ax.set_xlabel('G-test statistic')
    ax.set_ylabel('Distribution')
    ax.text(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(fams_lst)),
                horizontalalignment='right', verticalalignment='top', fontsize='large',
                transform=fig.transFigure)
    ax.set_title('Test cell type percentages '+oc+" gen2")
    sns.despine(fig=fig)
    fig.savefig(os.path.join(path_plot,'Test cell type percentages '+oc+' gen2.pdf'))


# In[25]:
    
#FIG 3D: frequency of cell type per generation 3

with open(os.path.join(path_pickle,'Fig3D'), 'rb') as fp:
    dct_plot = pickle.load(fp)
with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, cols_dct, xlabel='Progenitor', ylabel='Cells (%)', legend_title='Cell type',
                       legend_classes=class_dct['HSC'], legend_cols=cols_dct['HSC'])
    ax.set_title('Cell types percentages Generation 3')
    fig.savefig(os.path.join(path_plot,'Cell type percentages (Gen3).pdf'), bbox_inches='tight')
    
print('TEST GEN 3 CLASSES DEPENDENCE ON CONDITION')
def spread_support(support, data):
    res = np.zeros(len(support))
    d_sup,freq = np.unique(data, return_counts=True)
    res[[k in d_sup for k in support]] = freq
    return res

mc_iter = 250000.
for n_oc,oc in enumerate(or_cells):
    df_cnd = df_temp[(df_temp.Original_cell==oc)&(df_temp.Generation==3)]
    df_ary = np.array(df_cnd.loc[:,['Class','Family','Culture_condition']])
    cls_lst = np.unique(df_ary[:,0])
    fams_lst = np.unique(df_ary[:,1])
    fam_cl = np.array([df_ary[df_ary[:,1]==fam,0] for fam in fams_lst])
    fam_cnd = np.array([df_ary[df_ary[:,1]==fam,2][0] for fam in fams_lst])
    fam_data = np.array(list(zip(fam_cl, fam_cnd)))
    support = np.unique(df_ary[:,0])
    table = np.array([
            spread_support(support, np.hstack(fam_data[fam_data[:,1]==cnd,0]))
            for cnd in conds])
    stat = sps.chi2_contingency(table, correction=False, lambda_='log-likelihood')[0]
    per_indx = [np.random.permutation(fam_data[:,1]) for k in np.arange(mc_iter)]
    per_tables = np.array([
                [spread_support(support, np.hstack(fam_data[p_ind==cnd,0])) for cnd in conds]
            for p_ind in per_indx])
    per_stat = np.array([sps.chi2_contingency(p_tab, correction=False, lambda_='log-likelihood')[0]
            for p_tab in per_tables])
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    fig,ax=plt.subplots()
    sns.distplot(per_stat, norm_hist=True, kde=False, ax=ax)
    ax.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    ax.set_xlabel('G-test statistic')
    ax.set_ylabel('Distribution')
    ax.text(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(fams_lst)),
                horizontalalignment='right', verticalalignment='top', fontsize='large',
                transform=fig.transFigure)
    ax.set_title('Test cell type percentages '+oc+" gen3")
    sns.despine(fig=fig)
    fig.savefig(os.path.join(path_plot,'Test cell type percentages '+oc+' gen3.pdf'))


# In[26]:
    
#FIG 3E: frequency of cell type per generation 4

with open(os.path.join(path_pickle,'Fig3E'), 'rb') as fp:
    dct_plot = pickle.load(fp)
with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, cols_dct, xlabel='Progenitor', ylabel='Cells (%)', legend_title='Cell type',
                       legend_classes=class_dct['HSC'], legend_cols=cols_dct['HSC'])
    ax.set_title('Cell types percentages Generation 4')
    fig.savefig(os.path.join(path_plot,'Cell type percentages (Gen4).pdf'), bbox_inches='tight')

print('TEST GEN 4 CLASSES DEPENDENCE ON CONDITION')
def spread_support(support, data):
    res = np.zeros(len(support))
    d_sup,freq = np.unique(data, return_counts=True)
    res[[k in d_sup for k in support]] = freq
    return res

mc_iter = 250000.
for n_oc,oc in enumerate(or_cells):
    df_cnd = df_temp[(df_temp.Original_cell==oc)&(df_temp.Generation==4)]
    df_ary = np.array(df_cnd.loc[:,['Class','Family','Culture_condition']])
    cls_lst = np.unique(df_ary[:,0])
    fams_lst = np.unique(df_ary[:,1])
    fam_cl = np.array([df_ary[df_ary[:,1]==fam,0] for fam in fams_lst])
    fam_cnd = np.array([df_ary[df_ary[:,1]==fam,2][0] for fam in fams_lst])
    fam_data = np.array(list(zip(fam_cl, fam_cnd)))
    support = np.unique(df_ary[:,0])
    table = np.array([
            spread_support(support, np.hstack(fam_data[fam_data[:,1]==cnd,0]))
            for cnd in conds])
    stat = sps.chi2_contingency(table, correction=False, lambda_='log-likelihood')[0]
    per_indx = [np.random.permutation(fam_data[:,1]) for k in np.arange(mc_iter)]
    per_tables = np.array([
                [spread_support(support, np.hstack(fam_data[p_ind==cnd,0])) for cnd in conds]
            for p_ind in per_indx])
    per_stat = np.array([sps.chi2_contingency(p_tab, correction=False, lambda_='log-likelihood')[0]
            for p_tab in per_tables])
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    fig,ax=plt.subplots()
    sns.distplot(per_stat, norm_hist=True, kde=False, ax=ax)
    ax.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    ax.set_xlabel('G-test statistic')
    ax.set_ylabel('Distribution')
    ax.text(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(fams_lst)),
                horizontalalignment='right', verticalalignment='top', fontsize='large',
                transform=fig.transFigure)
    ax.set_title('Test cell type percentages '+oc+" gen4")
    sns.despine(fig=fig)
    fig.savefig(os.path.join(path_plot,'Test cell type percentages '+oc+' gen4.pdf'))


# In[27]:
    
#FIG 3F: frequency of cell type per generation 5

with open(os.path.join(path_pickle,'Fig3F'), 'rb') as fp:
    dct_plot = pickle.load(fp)
with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, cols_dct, xlabel='Progenitor', ylabel='Cells (%)', legend_title='Cell type',
                       legend_classes=class_dct['HSC'], legend_cols=cols_dct['HSC'])
    ax.set_title('Cell types percentages Generation 5')
    fig.savefig(os.path.join(path_plot,'Cell type percentages (Gen5).pdf'), bbox_inches='tight')

print('TEST GEN 5 CLASSES DEPENDENCE ON CONDITION')
def spread_support(support, data):
    res = np.zeros(len(support))
    d_sup,freq = np.unique(data, return_counts=True)
    res[[k in d_sup for k in support]] = freq
    return res

mc_iter = 250000.
for n_oc,oc in enumerate(or_cells):
    df_cnd = df_temp[(df_temp.Original_cell==oc)&(df_temp.Generation==5)]
    df_ary = np.array(df_cnd.loc[:,['Class','Family','Culture_condition']])
    cls_lst = np.unique(df_ary[:,0])
    fams_lst = np.unique(df_ary[:,1])
    fam_cl = np.array([df_ary[df_ary[:,1]==fam,0] for fam in fams_lst])
    fam_cnd = np.array([df_ary[df_ary[:,1]==fam,2][0] for fam in fams_lst])
    fam_data = np.array(list(zip(fam_cl, fam_cnd)))
    support = np.unique(df_ary[:,0])
    table = np.array([
            spread_support(support, np.hstack(fam_data[fam_data[:,1]==cnd,0]))
            for cnd in conds])
    stat = sps.chi2_contingency(table, correction=False, lambda_='log-likelihood')[0]
    per_indx = [np.random.permutation(fam_data[:,1]) for k in np.arange(mc_iter)]
    per_tables = np.array([
                [spread_support(support, np.hstack(fam_data[p_ind==cnd,0])) for cnd in conds]
            for p_ind in per_indx])
    per_stat = np.array([sps.chi2_contingency(p_tab, correction=False, lambda_='log-likelihood')[0]
            for p_tab in per_tables])
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    fig,ax=plt.subplots()
    sns.distplot(per_stat, norm_hist=True, kde=False, ax=ax)
    ax.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    ax.set_xlabel('G-test statistic')
    ax.set_ylabel('Distribution')
    ax.text(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(fams_lst)),
                horizontalalignment='right', verticalalignment='top', fontsize='large',
                transform=fig.transFigure)
    ax.set_title('Test cell type percentages '+oc+" gen5")
    sns.despine(fig=fig)
    fig.savefig(os.path.join(path_plot,'Test cell type percentages '+oc+' gen5.pdf'))
    
    
# In[28]:
    
#FIG 3G: frequency of cell type per generation 6

with open(os.path.join(path_pickle,'Fig3G'), 'rb') as fp:
    dct_plot = pickle.load(fp)
with sns.axes_style("ticks"):
    fig,ax = plot_bars(dct_plot, cols_dct, xlabel='Progenitor', ylabel='Cells (%)', legend_title='Cell type',
                       legend_classes=class_dct['HSC'], legend_cols=cols_dct['HSC'])
    ax.set_title('Cell types percentages Generation 6')
    fig.savefig(os.path.join(path_plot,'Cell type percentages (Gen5).pdf'), bbox_inches='tight')

print('TEST GEN 6 CLASSES DEPENDENCE ON CONDITION')
def spread_support(support, data):
    res = np.zeros(len(support))
    d_sup,freq = np.unique(data, return_counts=True)
    res[[k in d_sup for k in support]] = freq
    return res

mc_iter = 250000.
for n_oc,oc in enumerate(or_cells):
    df_cnd = df_temp[(df_temp.Original_cell==oc)&(df_temp.Generation==5)]
    df_ary = np.array(df_cnd.loc[:,['Class','Family','Culture_condition']])
    cls_lst = np.unique(df_ary[:,0])
    fams_lst = np.unique(df_ary[:,1])
    fam_cl = np.array([df_ary[df_ary[:,1]==fam,0] for fam in fams_lst])
    fam_cnd = np.array([df_ary[df_ary[:,1]==fam,2][0] for fam in fams_lst])
    fam_data = np.array(list(zip(fam_cl, fam_cnd)))
    support = np.unique(df_ary[:,0])
    table = np.array([
            spread_support(support, np.hstack(fam_data[fam_data[:,1]==cnd,0]))
            for cnd in conds])
    stat = sps.chi2_contingency(table, correction=False, lambda_='log-likelihood')[0]
    per_indx = [np.random.permutation(fam_data[:,1]) for k in np.arange(mc_iter)]
    per_tables = np.array([
                [spread_support(support, np.hstack(fam_data[p_ind==cnd,0])) for cnd in conds]
            for p_ind in per_indx])
    per_stat = np.array([sps.chi2_contingency(p_tab, correction=False, lambda_='log-likelihood')[0]
            for p_tab in per_tables])
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'wb') as fp:
        pickle.dump((per_stat,stat), fp)
    with open('./pickled_data/Cnd_diff_cond_uncorr_pop_gen2_'+oc, 'rb') as fp:
        per_stat,stat = pickle.load(fp)
    pval = ((per_stat>=stat).sum()+1.)/(mc_iter+1.)
    fig,ax=plt.subplots()
    sns.distplot(per_stat, norm_hist=True, kde=False, ax=ax)
    ax.vlines(stat, ymin=0., ymax=plt.gca().get_ylim()[1], colors='Red')
    ax.set_xlabel('G-test statistic')
    ax.set_ylabel('Distribution')
    ax.text(0.875, 0.875, 'p='+str(round(pval,6))+'\n'+'N. data='+str(len(fams_lst)),
                horizontalalignment='right', verticalalignment='top', fontsize='large',
                transform=fig.transFigure)
    ax.set_title('Test cell type percentages '+oc+" gen6")
    sns.despine(fig=fig)
    fig.savefig(os.path.join(path_plot,'Test cell type percentages '+oc+' gen6.pdf'))    


# In[29]:

# GENERATE STATISTICS FOR SUPPLEMENTARY TABLES

#SUPPLEMENTARY TABLE 1
for oc in or_cells:
    for cnd in conds:
        print(oc, cnd)
        df_temp = df[(df.Original_cell==oc)&(df.Culture_condition==cnd)]
        print('Total clones:', len(np.unique(df_temp.Family)))
        print('Clones in generation 0 clones:', len(np.unique(df_temp[df_temp.Generation==0].Family)))
        print('Clones with 2 cells in generation 1:',
              len([fam for fam in np.unique(df_temp[df_temp.Generation==1].Family)
                   if (df_temp[df_temp.Generation==1].Family==fam).sum()==2]))
        print('Clones recovered at 72h:', len(np.unique(df_temp[df_temp.Culture_time=='72h'].Family)))
        print('Cells in generation 1:', len(df_temp[df_temp.Generation==1]))
        print('Cells in generation 2:', len(df_temp[df_temp.Generation==2]))
        print('Cells in generation 3:', len(df_temp[df_temp.Generation==3]))
        print('Cells in generation 4:', len(df_temp[df_temp.Generation==4]))
        print('\n')
#SUPPLEMENTARY TABLE 3
for oc in or_cells:
    print(oc)
    for g in np.arange(4):
        df_temp = df[(df.Original_cell==oc)]
        print('Total cells in generation {0!s}: {1!s}'.format(g, (df_temp.Generation==g).sum()))
    print('Total cells in generation 4+: {!s}\n'.format((df_temp.Generation>=4).sum()))



