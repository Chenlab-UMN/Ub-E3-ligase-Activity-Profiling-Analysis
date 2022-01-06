# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:48:07 2021

@author: gy199
"""
import pandas as pd
import random
import statistics
import scipy.stats
import numpy as np
import os

        
#format  
#From	To
#Q92597	NDRG1
def read_uniprot2genename(file_source):
    file = open(file_source,"r")
    #line = file.readline()
    #headers = line.split('\t')
    linkdict = {}
    line = file.readline().strip('\n')
    while(line):
        contents = line.split(',') 
        linkdict[contents[0]] = contents[1]
        line = file.readline().strip('\n')   
    return(linkdict)

def read_E3dict(file_source):
    file = open(file_source,"r")
    line = file.readline()
    line = line.split(',')#skip header
    linkdict = {}
    line = file.readline().strip('\n')
    while(line):
        contents = line.split(',')
        subs = contents[1].split(';')
        linkdict[contents[0]] = subs
        line = file.readline().strip('\n')   
    return(linkdict)

def read_protein_groups(file_source):    
    proratiodict = {}
    proratiodf = pd.read_csv(file_source, delimiter = ",",header = 0)
    proratiodf = proratiodf[pd.to_numeric(proratiodf['ratio'], errors='coerce').notnull()]
    proratiodf.dropna(subset = ["ratio"], inplace=True)
    for i in proratiodf.index: 
            proratiodict[proratiodf.loc[i,'id']] = float(proratiodf.loc[i,'ratio'])
    return(proratiodict)


def read_ubiquitin_sites(file_source,proratiodict,unitogenedict,input_type,normbypro = False,log2trans = True):
    df = pd.read_csv(file_source, delimiter = ",",header = 0)
    df = df[pd.to_numeric(df['siteratio'], errors='coerce').notnull()]
    df.dropna(subset = ["siteratio"], inplace=True)
    
    if (normbypro==True):
        if(log2trans == True):
            for i in df.index:
                try:
                    proratio = proratiodict[df.loc[i,'id']]
                    normratio = float(df.loc[i,'siteratio'])/proratio
                    df.loc[i,'siteratio']= np.log2(normratio)
                except:
                    df.loc[i,'siteratio']= np.nan 
            
        else:
            for i in df.index:            
                try:
                    proratio = proratiodict[df.loc[i,'id']]
                    normratio = float(df.loc[i,'siteratio'])-proratio
                except:
                    df.loc[i,'siteratio']= np.nan 
            
                
    else:
        if(log2trans == True):
            for i in df.index:
                df.loc[i,'siteratio']= np.log2(df.loc[i,'siteratio'])
        else:
            pass
        
    df.dropna(subset = ["siteratio"], inplace=True)
    siteratiodf =  pd.DataFrame(columns = ['genename','proteinid','position','siteratio'])
    if (input_type=="p"):
        siteratiodf[['proteinid','position','siteratio']] = df[['id','position','siteratio']]
        for i in siteratiodf.index:    
            try:
                siteratiodf.loc[i,'genename']  = unitogenedict[siteratiodf.loc[i,'proteinid']]
            except:
                pass
    elif(input_type=="g"):
        siteratiodf[['genename','position','siteratio']] = df[['id','position','siteratio']]
    else:
        return()
    
    return(siteratiodf)


def sitetoavgproratio(df,input_type):
    avgprodf = pd.DataFrame(columns=['genename','proteinid','positions','siteratio','siteratios'])
    avgprorownum = 0
    if (input_type == "p"):
        for i in df.index: 
            idx = avgprodf.index[avgprodf['proteinid'] == df.loc[i,'proteinid']]
            if (len(idx) >0):
                    idx = idx[0]
                    avgprodf.loc[idx, 'positions'].append(df.loc[i,'position'])
                    if(df.loc[i,'siteratio']!=None):
                        avgprodf.loc[idx, 'siteratios'].append(df.loc[i,'siteratio'])
    
            elif(df.loc[i,'proteinid']!=''):
                    avgprodf.loc[avgprorownum,'genename']=df.loc[i,'genename']
                    avgprodf.loc[avgprorownum,'proteinid']=df.loc[i,'proteinid']
                    avgprodf.loc[avgprorownum,'positions']=[df.loc[i,'position']]
                    avgprodf.loc[avgprorownum,'siteratios']=[df.loc[i,'siteratio']]
                    avgprorownum = avgprorownum+1
                
    else:
        for i in df.index:
            idx = avgprodf.index[avgprodf['genename'] == df.loc[i,'genename']]
            if (len(idx) >0):
                    idx = idx[0]
                    avgprodf.loc[idx, 'positions'].append(df.loc[i,'position'])
                    if(df.loc[i,'siteratio']!=None):
                        avgprodf.loc[idx, 'siteratios'].append(df.loc[i,'siteratio'])
            elif(df.loc[i,'genename']!='NaN'):
                    avgprodf.loc[avgprorownum,'genename']=df.loc[i,'genename']
                    avgprodf.loc[avgprorownum,'proteinid']=df.loc[i,'proteinid']
                    avgprodf.loc[avgprorownum,'positions']=[df.loc[i,'position']]
                    avgprodf.loc[avgprorownum,'siteratios']=[df.loc[i,'siteratio']]
                    avgprorownum = avgprorownum+1
    avgprodf['siteratio'] = pd.DataFrame(avgprodf['siteratios'].values.tolist()).mean(1) 
    
    return(avgprodf)


def E3ligase_enrichment(siteratiodf,E3dict,prolevel):
    E3enrichdf =  pd.DataFrame(columns = ['e3ligase','site_count','site_avg','sample_avg','sample_std','p_value','substrates','substrate_ratio'])
    idx = 0
    ratio_sample = siteratiodf['siteratio']
    try:
        ratio_sample = [x for x in ratio_sample if ~np.isnan(x)]
    except:
        ratio_sample = list(filter(None, ratio_sample))

    for E3ligase in E3dict.keys():
        E3substrate = E3dict[E3ligase]
        try:
            subdf = siteratiodf.loc[siteratiodf['genename'].isin(E3substrate)]
        except:
            continue
        sitecount = len(subdf)
        substrate_list = []
        substrate_ratio_list = []
        if(sitecount<1):  #cannot calculate is there's no data points
            continue
        siteavg = statistics.mean(subdf['siteratio'])
        sampleavglist = []
        for i in range(0,10000):
            randsample = random.choices(ratio_sample, k=sitecount)
            sampleavglist.append(statistics.mean(randsample))
        sampleavg = statistics.mean(sampleavglist)
        samplestd = statistics.stdev(sampleavglist)
        zscore = abs(siteavg-sampleavg)/samplestd
        pvalue = scipy.stats.norm.sf(abs(zscore))*2
        
        
        if(prolevel == True):
            for i in subdf.index:
                if (subdf.loc[i,'genename'] not in substrate_list):
                     substrate_list.append(subdf.loc[i,'genename'])
                     substrate_ratio_list.append(subdf.loc[i,'genename']+':'+str(np.round(subdf.loc[i,'siteratio'],5)))

        else:
            for i in subdf.index:
                if (subdf.loc[i,'genename'] not in substrate_list):
                     substrate_list.append(subdf.loc[i,'genename'])
                try:
                    substrate_ratio_str = subdf.loc[i,'genename']+'_'+str(int(subdf.loc[i,'position']))
                except:
                    substrate_ratio_str = subdf.loc[i,'genename']+'_'+"NaN"
                if (substrate_ratio_str not in substrate_ratio_list):
                     substrate_ratio_list.append(substrate_ratio_str+':'+str(np.round(subdf.loc[i,'siteratio'],5)))
 
    
        E3enrichdf.loc[idx] = [E3ligase,sitecount,siteavg,sampleavg,samplestd,pvalue,substrate_list,substrate_ratio_list]
        idx = idx+1
    return(E3enrichdf)



def group_siteratiodf(df):
    del_index = [] 
    df["leading_e3ligase"] =""
    df["e3ligase_group"] =""
    for i in df.index:
        df.at[i,"e3ligase_group"] = [df.loc[i,"e3ligase"]]
        df.at[i,"leading_e3ligase"] = [df.loc[i,"e3ligase"]]

    for i in df.index:
        if (i not in del_index):
            for j in df.index:
                if(i!=j and j not in del_index):
                    if (set(df.loc[j,"substrates"]).issubset(set(df.loc[i,"substrates"]))):  
                        if (set(df.loc[i,"substrates"]).issubset(set(df.loc[j,"substrates"]))):
                            df.loc[i,"leading_e3ligase"].extend(df.loc[j,"leading_e3ligase"])
                            df.loc[i,"e3ligase_group"].extend(df.loc[j,"e3ligase_group"])
                            del_index.append(j)
                        else:
                            df.loc[i,"e3ligase_group"].extend(df.loc[j,"e3ligase_group"])
                            del_index.append(j)                    
    df = df.drop(del_index)
    df = df.drop(columns=['e3ligase'])
    df = df[["leading_e3ligase", "e3ligase_group", "site_count", "site_avg", "sample_avg",	"sample_std", "p_value", "substrates",	"substrate_ratio"]]
    return(df)    
    
    
def e3enrich(siteratio_dir,input_type,output_dir,exp_label='',proratio_dir="None",grouped=False,ratio_output=False,log2trans=True):     
    #input_type: 'UniprotAC' or "genename"
    module_dir = os.path.abspath(os.path.dirname(__file__))
    unitogenedict = read_uniprot2genename(module_dir+'/unigenedict.csv') 
    E3dict = read_E3dict(module_dir+'/E3dict.csv') 
    
    
    #############################################################################################
    #input_type section, return() if input_type is invalid
    if (input_type=="UniprotAC" or input_type=="Uniprot" or input_type=="protein"):
        input_type = "p"
    elif (input_type=="gene symbol" or input_type=="gene name" or input_type=="gene"):
        input_type = 'g'
    else:
        print("Invalid input_type, please use one of the following strings")
        print("\"UniprotAC\"or \"protein\"  example:Q00987,P40337,Q9HAU4,O43791,Q7Z6Z7")
        print("\"gene symbol\"or \"gene\"  example:MDM2,VHL,SMUF2,SPOP,HUWE1")
        return()
   
    
    print("Analyzing data of experiment "+exp_label)
    #############################################################################################
    #gather ratios from input filesand perform enrichment analysis   
    if(proratio_dir=="None"):
        normbypro = False
        proratiodict = {}
    else:
        normbypro = True
        proratiodict = read_protein_groups(proratio_dir)
    siteratiodf = read_ubiquitin_sites(siteratio_dir,proratiodict,unitogenedict,input_type,normbypro,log2trans)
    avgsiteratiodf = sitetoavgproratio(siteratiodf,input_type) #from site ratio to average protein ratio
    
    proE3enrichdf = E3ligase_enrichment(avgsiteratiodf,E3dict,True)
    siteE3enrichdf = E3ligase_enrichment(siteratiodf,E3dict,False)
    
    #############################################################################################
    #export site level ratio file, protein level ratio file and E3_ubiprotein linkage file  
    if(ratio_output==True):
        siteratiodf.to_csv(output_dir+'/site_level_ratio'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8') 
        outavgsiteratiodf = avgsiteratiodf
        for i in outavgsiteratiodf.index: 
            outavgsiteratiodf.loc[i, 'number_of_sites'] = len(outavgsiteratiodf.loc[i,'positions'])
            string1 = []
            for float in outavgsiteratiodf.loc[i, 'positions']:
                if(np.isnan(float)):
                    string1.append("NaN")
                else:
                    string1.append(str(int(float)))
            #string1 = [str(int(float)) for float in outavgsiteratiodf.loc[i, 'positions']]
            outavgsiteratiodf.loc[i,'positions'] = ';'.join(string1)
            string2 = [str(float) for float in outavgsiteratiodf.loc[i, 'siteratios']]
            outavgsiteratiodf.loc[i,'siteratios'] = ';'.join(string2)
            outavgsiteratiodf = outavgsiteratiodf.rename(columns={"siteratio": "average_siteratio"})
        outavgsiteratiodf.to_csv(output_dir+'/protein_level_ratio'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8') 
        
        linkdf =  pd.DataFrame(columns = ['e3ligase','substrate_protein'])
        linkidx = 0
        for i in proE3enrichdf.index:
            for j in proE3enrichdf.loc[i,'substrates']:
                linkdf.loc[linkidx,'e3ligase'] = proE3enrichdf.loc[i,'e3ligase']
                linkdf.loc[linkidx,'substrate_protein'] = j
                linkidx = linkidx+1
        linkdf.to_csv(output_dir+'/e3_sub_network'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8')    
        
    #############################################################################################
    #organize E3ligase enrichment results based on different parameters 
    #para1: grouped or not
    #para2: norm by protein level or not


    if (grouped == True):
        proE3enrichdf = group_siteratiodf(proE3enrichdf)
        siteE3enrichdf = group_siteratiodf(siteE3enrichdf)
        for i in proE3enrichdf.index:  
            proE3enrichdf.loc[i,'leading_e3ligase'] = ';'.join(proE3enrichdf.loc[i,'leading_e3ligase'])
            proE3enrichdf.loc[i,'e3ligase_group'] = ';'.join(proE3enrichdf.loc[i,'e3ligase_group'])
            proE3enrichdf.loc[i,'substrates'] = ';'.join(proE3enrichdf.loc[i,'substrates'])
            proE3enrichdf.loc[i,'substrate_ratio'] = ';'.join(proE3enrichdf.loc[i,'substrate_ratio'])

        for i in siteE3enrichdf.index: 
            siteE3enrichdf.loc[i,'leading_e3ligase'] = ';'.join(siteE3enrichdf.loc[i,'leading_e3ligase'])
            siteE3enrichdf.loc[i,'e3ligase_group'] = ';'.join(siteE3enrichdf.loc[i,'e3ligase_group'])
            siteE3enrichdf.loc[i,'substrates'] = ';'.join(siteE3enrichdf.loc[i,'substrates'])
            siteE3enrichdf.loc[i,'substrate_ratio'] = ';'.join(siteE3enrichdf.loc[i,'substrate_ratio'])

        if (normbypro == True):
                proE3enrichdf.to_csv(output_dir+'/E3enrichment_protein_level_normalized_by_protein_grouped'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8') 
                siteE3enrichdf.to_csv(output_dir+'/E3enrichment_site_level_normalized_by_protein_grouped'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8') 
        else:   
                proE3enrichdf.to_csv(output_dir+'/E3enrichment_protein_level_grouped'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8') 
                siteE3enrichdf.to_csv(output_dir+'/E3enrichment_site_level_grouped'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8')    
    else:
        for i in proE3enrichdf.index:  
            proE3enrichdf.loc[i,'substrates'] = ';'.join(proE3enrichdf.loc[i,'substrates'])
            proE3enrichdf.loc[i,'substrate_ratio'] = ';'.join(proE3enrichdf.loc[i,'substrate_ratio'])

        for i in siteE3enrichdf.index:  
            siteE3enrichdf.loc[i,'substrates'] = ';'.join(siteE3enrichdf.loc[i,'substrates'])
            siteE3enrichdf.loc[i,'substrate_ratio'] = ';'.join(siteE3enrichdf.loc[i,'substrate_ratio'])

        if (normbypro == True):
                proE3enrichdf.to_csv(output_dir+'/E3enrichment_protein_level_normalized_by_protein'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8') 
                siteE3enrichdf.to_csv(output_dir+'/E3enrichment_site_level_normalized_by_protein'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8') 
        else:   
                proE3enrichdf.to_csv(output_dir+'/E3enrichment_protein_level'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8') 
                siteE3enrichdf.to_csv(output_dir+'/E3enrichment_site_level'+str(exp_label)+'.csv', sep=',',index = False, encoding='utf-8')  
    return()



 




 