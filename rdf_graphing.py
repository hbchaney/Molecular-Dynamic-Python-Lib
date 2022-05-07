# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 17:00:27 2021

@author: Harrison
"""

import matplotlib.pyplot as plt 


def file_open(filepath):
    '''
    takes a files complete name with path included 
    '''
    
    with open(filepath) as f: 
        lines = f.readlines()
        
    # print(lines)
    data1 = lines[4:504]
    data1 = [data.split(' ') for data in data1]
    print(data1)
    
    x = []
    y = []
    for values in data1:
        x.append(float(values[1]))
        y.append(float(values[2]))
        
    plt.plot(x,y,'.-')
    plt.show()
        
    
file_open('G:/Shared drives/1 Kathy Lu Group/Harrison/reaxff/midcarbon/rdf/t7_shd.rdf')

