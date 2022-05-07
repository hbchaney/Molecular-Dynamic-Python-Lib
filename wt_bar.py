# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:04:42 2022

@author: hbcha
"""


import json 
import matplotlib.pyplot as plt
import numpy as np 

def json_wt_bar (filename,polymer_name) : 
    
    with open(filename,'r') as f: 
        con = json.load(f)
    
    start = con[0] 
    end = con[-1] 
    
    print(start['other']['num'])
    
    Si1,O1,C1,H1 = start['other']['num']
    Si2,O2,C2,H2 = end['other']['num']
    
    
    
    total1 = Si1*28.085 + O1*16 + C1*12.01 + H1 
    total2 = Si2*28.085 + O2*16 + C2*12.01 + H2
    
    
    Si = np.array([Si1*28.085/total1,Si2*28.085/total2])
    O = np.array([O1*16/total1,O2*16/total2])
    C = np.array([C1*12.01/total1,C2*12.01/total2])
    H = np.array([H1/total1,H2/total2])
    
    # print(Si,O,C,H)
    
    width = .35
    
    labels = [f'{polymer_name} start',f'{polymer_name} end']
    
    fig, ax = plt.subplots()
    
    location = [1,1.7]
    
    ax.bar(location,Si,width,tick_label=labels,label='Si',align='center')
    ax.bar(location,O,width,tick_label=labels,bottom=Si,label='O',align='center')
    ax.bar(location,C,width,tick_label=labels,bottom=Si+O,label='C',align='center')
    ax.bar(location,H,width,tick_label=labels,bottom=Si+O+C,label='H',align='center')
    
    ax.set_ylabel('weight fraction')
    ax.set_title(f'changes in compositions for {polymer_name}')
    ax.legend(loc='center')
    
    plt.show() 
    
    
json_wt_bar(r'G:\My Drive\Research\reaxff\PSO\gas_results\mc-t11-gas.json','PSO')
    