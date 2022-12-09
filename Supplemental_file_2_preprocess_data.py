##########################################################################################
### Overview: This script is to load the data and packages that you need for the analysis
###           It will also load in the experimental data and process and merge for downstream
###           plotting and statistical analysis
###
### Correspondence: Leila Perie (leila.perie@curie.fr) ; Alessandro Donada alessandro.donada@gmail.com)
###
### Date: 09/12/2022
###########################################################################################


#!/usr/bin/env python
# coding: utf-8

# In[1]: 0_process_data
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
path_plot = os.path.join(path_proj, 'figures')


# In[2]:

#(USER-DEFINED)     

#when we make our  plots we want each  cell-type to plotted using a different color. here  we define a dictionary called cell_cols that returns a different color for each celltype.
cell_cols = {'HSC':'seagreen', 'MPP':'Limegreen', 'LMPP':'Red', 'CMP':'Deepskyblue','GMP':'Violet', 'MEP':'Royalblue', 'CD34-':'aqua'}
sns.palplot(cell_cols.values())


#  cell_class_exp_time is a function that determines a cell's class upon its markers.
#  description: in the multigen experiments we use protein expression on the cell surface to define
#               a cell's identity. This function takes as input a dataframe with  protein  expression for each cell
#               and threshold values for each protein. If a cells protein measurement exceeds some threshold
#               value we say that the cell is positive for that protein. Based on which  proteins the cells is
#                expressing we can define what cell type it should be.
#
# inputs: df = protein expression dataframe
# inputs: thr_ = a threshold value for the different proteins
# outputs: celltype, color, rank. The rank is used for plotting ordering

def cell_class_exp_time(df, thr_cd34, thr_cd38, thr_cd90, thr_cd45l, thr_cd45h, thr_cd123):

    if df['CD34'] > thr_cd34:
        if df['CD38'] <= thr_cd38:
            if df['CD90'] > thr_cd90:
                return {'class':'HSC', 'color':cell_cols['HSC'], 'rank':1}
            else:
                if df['CD45'] > thr_cd45h:
                    return {'class':'LMPP', 'color':cell_cols['LMPP'], 'rank':2}
                else:
                    return {'class':'MPP', 'color':cell_cols['MPP'], 'rank':3}
        else:
            if df['CD123'] > thr_cd123:
                if df['CD45'] > thr_cd45l:
                    return {'class':'GMP', 'color':cell_cols['GMP'], 'rank':4}
                else:
                    return {'class':'CMP', 'color':cell_cols['CMP'], 'rank':5}
            else:
                return {'class':'MEP', 'color':cell_cols['MEP'], 'rank':6}
    
    else:
        return {'class':'CD34-', 'color':cell_cols['CD34-'], 'rank':7}

    return {'class':'NC', 'color':'Black', 'rank':100}#'NC'


#this function ensures that the dataframe contains information related to the Experiment name and the cell culture length.
def cell_class(df):
    if df.Experiment+'_'+df.Culture_time in dct_thr.keys():
        return cell_class_exp_time(df, *(dct_thr[df.Experiment+'_'+df.Culture_time]))
    else:
        return {'class':'Experiment_or_time_not_found', 'color':'Black', 'rank':100}#'Experiment_or_time_not_found'    

#this line imports the threshold values defined when analysing the bulk wells
df_gating = pd.read_excel(os.path.join(path_proj,'Gating_matrix.xlsx'), index_col=0)
dct_thr = {exp_time:df_gating[exp_time].values for exp_time in df_gating.columns}
dct_thr


# In[3]:

#(EXPERIMENT SPECIFIC)    
#conds and cond_rule are used later to rename the culture conditions


#This line add the experimental conditions tested, in particular the different cell culture cocktails used
conds = ['GT', 'Diff'] 
def cond_rule(cond):
    if cond == 'GT' :
        return conds[0]
    elif cond == 'Diff':
        return conds[1]
    else:
        return 'NA'


#this line set the function used to detect the impossible families, defined as families that display too many cells for the theoretical limit, which is established
#based on the number of divisions performed: a family with cells that divided 2 times cannot display more than 2^2 cells
def find_impossible_families(df):
    fams = np.unique(df.Family)
    lst_fam_vec = [np.unique(df[df.Family==fam].Generation, return_counts=True) for fam in fams]
    cohort_number = np.array([(el[1]/np.power(2.,el[0])).sum() for el in lst_fam_vec])
    if all(cohort_number<=1):
        print('No impossible families detected')
        return []
    else:
        impossible_families = fams[cohort_number>1]
        print('Impossible families:', impossible_families)
        return impossible_families


# In[4]:

#Import the single cell data from the flow cytometry experiment

lst_sc_files = [file for file in os.listdir(path_sc) if '.csv' in file and file!='Pooled_data.csv']
df_sc_lst = []


# iterate over all of the different files merging and processing them.
# the processed and pooled data is saved to Pooled_data.csv
for file_name in lst_sc_files:
    print(file_name)
    
#(EXPERIMENT SPECIFIC) as the filenames here ends with '.csv'
    experiment_name = file_name[:-4]
    
#(EXPERIMENT SPECIFIC) sep=',', decimal='.' are unusual
    df = pd.read_csv(os.path.join(path_sc, file_name), sep=',', decimal='.')
    print(df.columns)
    
#Here the user should decide how to format the imported data
    df = df.assign(
        Culture_condition=df.apply(lambda r: cond_rule(r.Condition), axis=1),
        Experiment=[experiment_name for k in range(len(df))]
    )
    df = df.assign(
        Well_experiment=lambda r:r.Well+'_'+r.Experiment,
        Family=lambda r: r.Well + r.Color + r.Culture_time + r.Original_cell + r.Experiment,
        Class=df.apply(func=lambda r:cell_class(r)['class'], axis=1),
        Cell_color=df.apply(func=lambda r:cell_class(r)['color'], axis=1),
        Cell_rank=df.apply(func=lambda r:cell_class(r)['rank'], axis=1)
    )
   
#Detect and remove impossible families
    impossible_fams = find_impossible_families(df)
    df = df[~(df.Family.isin(impossible_fams))]
    df_sc_lst.append(df)
  
df_pool = pd.concat(df_sc_lst, ignore_index=True, sort=True)


#Export processed and pooled data
df_pool.to_csv(os.path.join(path_sc)+'\Pooled_data.csv', sep=',', decimal='.', index=False)


# In[5]:
   

np.unique(df_pool.Class, return_counts=True), np.unique(df_pool.Culture_time, return_counts=True)


