#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:59:09 2017

@author: na399
"""

import pandas as pd
from sklearn.preprocessing import scale
import seaborn as sns

sns.set_style("whitegrid")



df_CAGE = pd.read_table("htbinf_cage_tpms.txt")

tissue_names = df_CAGE.columns[3:10]

array_CAGE_scaled = scale(df_CAGE.iloc[:, 3:10], axis=1)

df_CAGE_scaled = pd.DataFrame(array_CAGE_scaled, 
                             columns=tissue_names +'_scaled')

df_CAGE_2 = pd.concat([df_CAGE, df_CAGE_scaled], axis=1)

df_CAGE_plot = df_CAGE_2[['tc_id', 'cer', 'emb', 'liv', 'lun', 'mac', 'som', 'vis']]
df_CAGE_plot_scaled = df_CAGE_2[['tc_id',  'cer_scaled', 'emb_scaled', 'liv_scaled', 'lun_scaled', 'mac_scaled', 'som_scaled', 'vis_scaled']]

df_CAGE_plot_melted = pd.melt(df_CAGE_plot.iloc[0:30,], id_vars='tc_id', var_name='tissue', value_name='tpm')
df_CAGE_plot_scaled_melted = pd.melt(df_CAGE_plot_scaled.iloc[0:30,], id_vars='tc_id', var_name='tissue', value_name='tpm')


ax = sns.swarmplot(x="val", y="tc_id", data=df_CAGE_plot_melted)

ax = sns.swarmplot(x="val", y="tc_id", data=df_CAGE_plot_scaled_melted)