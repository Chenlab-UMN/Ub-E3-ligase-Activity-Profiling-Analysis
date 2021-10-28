# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 15:23:21 2021

@author: gong
"""

import ube3_apa


working_dir = "C:/Users/gy199/OneDrive/Desktop/E3_ligase_enrichment/"
ube3_apa.e3enrich(working_dir+"testdata/siteratio_testdata.csv","UniprotAC",working_dir+"testdata",exp_label="1",grouped=False,ratio_output=True,log2trans=False)
ube3_apa.e3enrich(working_dir+"testdata/siteratio_testdata.csv","UniprotAC",working_dir+"testdata",exp_label="2",grouped=True,ratio_output=True,log2trans=False)
ube3_apa.e3enrich(working_dir+"testdata/siteratio_testdata.csv","UniprotAC",working_dir+"testdata",exp_label="3",proratio_dir=working_dir+"testdata/proteinratio_testdata.csv",grouped=False,ratio_output=True,log2trans=False)



