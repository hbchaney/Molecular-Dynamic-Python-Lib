# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 11:11:02 2022

@author: hbcha
"""

import os 
import pandas as pd 
import re 

#test
# cwd = os.getcwd()

# os.chdir(r'G:\My Drive\Gold Nano Results')
# print(os.getcwd())

def file_parse(direct):
    """

    Parameters
    ----------
    direct : str 
        put in the directory where all the folders are stored 

    Returns
    -------
    a pandas data frame with the info set up 

    """
    
    #make sure your directory with all the polymer names has only that!!!!
    os.chdir(direct)
    polymers = os.listdir()
    # print(polymers)
    cwd = os.getcwd()
    
    polymer_d = {}
    
    for poly in polymers:
        os.chdir(poly)
        # print(os.getcwd())
        
        ply_df = pd.DataFrame()
        # print(type(ply_df))
        # ply_df.loc['2nm','0.01'] = 6
        # print(ply_df)
        
        for files in os.listdir(): 
            # print(files)
            nm, vol = re.findall(r'(\d+nm)\D*(\d+\.\d+).*%',files)[0]
            # print(nm,vol)
            ply_df.loc[nm,vol] = 'hello'
            
        
        for files in os.listdir():
            nm, vol = re.findall(r'(\d+nm)\D*(\d+\.\d+).*%',files)[0]
            # print(nm,vol)
            Data = YiJe(files)
            ply_df.at[nm,vol] = Data
            
            

                    
        polymer_d[poly] = ply_df
        
        
        
        os.chdir(cwd)
    
    return polymer_d
# file_parse(r'G:\My Drive\Gold Nano Results')
    

def YiJe(file_name):
    """
    written by Yi Je 

    Parameters
    ----------
    file_name : str
        name of the file

    Returns
    -------
    max

    """
    
    sheet = pd.read_csv(file_name, skiprows=4).iloc[0:100,1:10]
    sheet['mean'] = sheet.iloc[:,1:10].mean(axis=1)
    max_value = pd.DataFrame(sheet.loc[sheet['mean'].idxmax()]).T
    
    return max_value


# print(YiJe(r'G:\My Drive\Gold Nano Results\PVA\PVA_2nm_vol0.50%_Results.csv'))
    
x = file_parse(r'G:\My Drive\Gold Nano Results') # put in your directory with all the polymer info
print(x['PDMS'].loc[:,'10.0']) #shows the dframe data for pdms for every size at 10% vol 
# print(x['PS'].loc['2nm',:]) #shows the dframe data for ps at 2nm for every vol 
    
    
    
    
    
    