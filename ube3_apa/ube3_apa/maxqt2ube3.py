# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 15:01:52 2021

@author: gy199
"""

import pandas as pd
import pathlib

def read_modification(file_source,output_dir,exp_label = ""):    
    file = open(file_source,"r")
    line = file.readline()
    headers = line.split('\t')
    if (exp_label == ""):
        column_header = 'Ratio'
    else:
        column_header = 'Ratio '+exp_label
    
    df = pd.DataFrame(columns=['id','position','siteratio'])
    tempdict = {}
    idindex = 0   
    line = file.readline()
    while(line):
        contents = line.split('\t')
        for i in range(0,len(headers)):
            tempdict[headers[i]] = contents[i]
        if (tempdict["Reverse"]!= '+' and tempdict["Potential contaminant"]!='+' and float(tempdict["Localization prob"])>=0.75 and tempdict["Protein"]!=''):
                proid =  tempdict['Protein'].strip('\t').split('|')[1]
                df.loc[idindex] = [proid,tempdict['Position'],tempdict[column_header]] 
                idindex = idindex+1
        else:
            pass
        tempdict = {}
        line = file.readline()
    exp_label = exp_label.replace("/", "")
    path = pathlib.Path(output_dir,'site_'+str(exp_label)+'.csv')
    df.to_csv( path_or_buf=path, sep=',', index=False, mode='w')
    return(path)


def read_modification_normalized(file_source,output_dir,exp_label = ""):    
    file = open(file_source,"r")
    line = file.readline()
    headers = line.split('\t')
    if (exp_label == ""):
        column_header = 'Ratio normalized'
    else:
        column_header = 'Ratio '+exp_label+' normalized'
    
    df = pd.DataFrame(columns=['id','position','siteratio'])
    tempdict = {}
    idindex = 0   
    line = file.readline()
    while(line):
        contents = line.split('\t')
        for i in range(0,len(headers)):
            tempdict[headers[i]] = contents[i]
        if (tempdict["Reverse"]!= '+' and tempdict["Potential contaminant"]!='+' and float(tempdict["Localization prob"])>=0.75 and tempdict["Protein"]!=''):
                proid =  tempdict['Protein'].strip('\t').split('|')[1]
                df.loc[idindex] = [proid,tempdict['Position'],tempdict[column_header]] 
                idindex = idindex+1
        else:
            pass
        tempdict = {}
        line = file.readline()
    exp_label = exp_label.replace("/", "")
    path = pathlib.Path(output_dir,'sitenorm_'+str(exp_label)+'.csv')
    df.to_csv( path_or_buf=path, sep=',', index=False, mode='w')
    return(path)

def read_proteingroup(file_source,output_dir,exp_label = ""):    
    file = open(file_source,"r")
    line = file.readline()
    headers = line.split('\t')
    if (exp_label == ""):
        column_header = 'Ratio'
    else:
        column_header = 'Ratio '+exp_label
    
    df = pd.DataFrame(columns=['id','ratio'])
    tempdict = {}
    idindex = 0   
    line = file.readline()
    while(line):
        contents = line.split('\t')
        for i in range(0,len(headers)):
            tempdict[headers[i]] = contents[i]
        if (tempdict["Reverse"]!= '+' and tempdict["Potential contaminant"]!='+' ):
                proids = tempdict['Majority protein IDs'].strip('\t').split(';')
                for proid in proids:
                    try:
                        proid =  proid.split('|')[1]
                        df.loc[idindex] = [proid,tempdict[column_header]] 
                        idindex = idindex+1
                    except:
                        pass
        else:
            pass
        tempdict = {}
        line = file.readline()
    exp_label = exp_label.replace("/", "")
    path = pathlib.Path(output_dir,'pro_'+str(exp_label)+'.csv')
    df.to_csv( path_or_buf=path, sep=',', index=False, mode='w')
    return(path)

def maxqt2ube3(output_dir,site_ratio_dir,pro_ratio_dir='',exp_label=''):
    print('Converting MaxQuant results of Experiment '+str(exp_label)+' to standard input')
    if(pro_ratio_dir!=''):
        pathoriginal = read_modification(site_ratio_dir,output_dir,exp_label)
        pathpro = read_proteingroup(pro_ratio_dir,output_dir,exp_label)
        return(pathoriginal,pathpro)
    else:
        pathnorm = read_modification_normalized(site_ratio_dir,output_dir,exp_label)
        return(pathnorm,'None')
        
        



