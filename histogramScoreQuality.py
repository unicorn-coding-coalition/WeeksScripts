# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 00:02:55 2019
Create Histograms of Fpocket Score against whether it is a pocket or not
@author: nlama
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
os.chdir("D:/Weeks/PocketFinding/dataSetA")
###Plot Parameters
red_patch = mpatches.Patch(color='indianred', label='Not a Pocket')
blue_patch = mpatches.Patch(color='steelblue',label='Pocket')


inDf = pd.read_csv("defaultFpocketDescriptors.csv")

inDf.plot(kind='bar',x=['PDB','Pocket'],y='Score',color=[np.where(inDf["Quality"]==1, 'steelblue', 'indianred')])
plt.legend(handles=[blue_patch,red_patch])
inDf.sort_values(by='Score', ascending=False,inplace=True)
inDf.plot(kind='bar',x=['PDB','Pocket'],y='Score',color=[np.where(inDf["Quality"]==1, 'steelblue', 'indianred')])
plt.legend(handles=[blue_patch,red_patch])

#inDf.plot(subplots=True,kind='bar',x=['PDB','Pocket'],color=[np.where(inDf["Known Binding Site"]==1, 'steelblue', 'indianred')], layout=(4,4))
#inDf.plot(kind='bar',x=['PDB','Pocket'], layout=(4,4))

#
#for col in inDf.columns:
#    try:
#        plt.rcParams.update({'font.size': 22})
#        inDf.plot(kind='bar',x=['PDB','Pocket'],y=col,color=[np.where(inDf["Known Binding Site"]==1, 'steelblue', 'indianred')])
#        plt.title(col)
#        plt.legend(handles=[blue_patch,red_patch])
#        plt6.savefig('{0}Hist.pdf'.format(col))
#        print(col, inDf[col].corr(inDf['Quality'])) ##Print correlations of each column to quality
#    except:
#        continue


#plot.hist(inDf['Quality'],inDf['Score'])
