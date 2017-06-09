#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:59:09 2017

@author: na399
"""

import pandas as pd
from sklearn.preprocessing import scale

df_CAGE = pd.read_table("htbinf_cage_tpms.txt")

tissue_names = df_CAGE.columns[3:10]

array_CAGE_scaled = scale(df_CAGE.iloc[:, 3:10], axis=1)

df_CAGE_scaled = pd.DataFrame(array_CAGE_scaled, 
                             columns=tissue_names +'_scaled')

df_CAGE_2 = pd.concat([df_CAGE, df_CAGE_scaled], axis=1)

threshold = 0.0


def isAbovethreshold(value):
    return value >= threshold


for tissue in tissue_names:
    df_CAGE_2[tissue + '_threshold'] = map(isAbovethreshold, df_CAGE_2[tissue + '_scaled'])


def isTissueSpecific_threshold(tissue_threshold, other_tissues_threshold):
    if tissue_threshold == True:
        if sum(other_tissues_threshold) == 0:
            return True
        else:
            return False
    else:
        return False


for tissue in tissue_names:
    other_tissues = tissue_names.delete(tissue_names.tolist().index(tissue))
    for index, row in df_CAGE_2.iterrows():
        df_CAGE_2.loc[index, tissue + '_t_specific'] = isTissueSpecific_threshold(
            df_CAGE_2.loc[index, tissue + '_threshold'], df_CAGE_2.loc[index, other_tissues + '_threshold'])

for tissue in tissue_names:
    print df_CAGE_2[tissue + '_t_specific'].sum()

df_CAGE_3 = df_CAGE_2[['location', 'strand', 'cer_t_specific', 'emb_t_specific',
                       'liv_t_specific', 'lun_t_specific', 'mac_t_specific',
                       'som_t_specific', 'vis_t_specific']]

for index, row in df_CAGE_3.iterrows():
    if sum(row[-7:]) == 0:
        df_CAGE_3.drop(index, inplace=True)

df_CAGE_3['tissue'] = ''

for tissue in tissue_names:
    for index, row in df_CAGE_3.iterrows():
        if row[tissue + '_t_specific'] == 1:
            df_CAGE_3.loc[index, 'tissue'] = tissue

df_CAGE_4 = df_CAGE_3[['location', 'strand', 'tissue']]

df_CAGE_4.to_csv("cage_tissue.csv")
