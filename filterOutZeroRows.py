# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 16:29:56 2018
This script is made to read in an edge list and filter out any zeros.

Note: Should update this to reformat ct file so I don't have to do it 
manually again.

@author: nlama
"""

import pandas as pd

df = pd.read_csv("D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles/rnpb2-Copy.csv",
                 header = "infer")

df_filt = df.loc[df['Target'] != 0]

pd.DataFrame.to_csv(df_filt, "D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles/rnpb2Filt.csv",
                    columns= ["Source","Target"],index = False)