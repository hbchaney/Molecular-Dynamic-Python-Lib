# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 10:11:24 2021

@author: Harrison
"""


import os 
import matplotlib.pyplot as plt
# print(os.listdir())
# os.chdir('D:/College Work/Research/Polymer_Derived/python_gas/Si_t3')    #change dir here if bonds are in a different folder 
# print(os.listdir())


def Si_Simultidir(): 
    os.chdir('G:/My Drive/Research/reaxff/hc/hc-long-gas')
    dir_present = os.listdir()
    dir_filter = []
    for directs in dir_present: 
        # print(directs[-3])
        if directs[-3] != '.':
            dir_filter.append(directs)
    
    # print(dir_filter)
    cwd = os.getcwd()
    Si_Si = []
    Si_C = []
    
    for directs in dir_filter: 
        os.chdir(directs)
        bonds = Si_Sigroup()
        Si_Si = Si_Si + bonds[0]
        Si_C =  Si_C + bonds[1]
        os.chdir(cwd)
    
    ts = range(0,len(Si_Si))
    
    # fig, axs = plt.subplots(2,1,constrained_layout=True)
    # fig.suptitle('Si-Si vs Si-C')
    # axs[0].plot(ts,Si_Si)
    # axs[0].set_title('Si-Si bonds')
    # axs[1].plot(ts,Si_C)
    # axs[1].set_title('Si-C bonds')
    
    plt.plot(ts,Si_Si,label = 'Si-Si bonds')
    plt.plot(ts,Si_C,label = 'Si-C bonds')
    plt.legend()
    plt.show()
        
    

    
    

def Si_Sigroup(): 
    files = os.listdir()
    # print(sorted(files))
    # print(files[0][13:-6])
    files = sorted(files, key=lambda x: float(x[13:-6]))
    # print(files)
    
    Si_Si = []
    Si_C =[]
    
    for file in files: 
        bonds_Si_Si = Si_Sisingle(file)[0]
        bonds_Si_C = Si_Sisingle(file)[1]
        Si_Si.append(bonds_Si_Si)
        Si_C.append(bonds_Si_C)
        
    # print(Si_Si)
    return [Si_Si, Si_C]
       
    # fig, axs = plt.subplots(2)
    # fig.suptitle('Vertically stacked subplots')
    # axs[0].plot(ts,Si_Si)
    # axs[1].plot(ts,Si_C)
     
     
     
     


    
    
# 'gas_creation.1.reaxc'   
def Si_Sisingle(file_name): 
    '''
    takes in a bond info file and finds the Si bond info for the first time step in the file
    outputs a list of Si and Si bonds 
    
    the number of Si_Si bonds found 
    
    
    '''
    
    with open(file_name) as f: 
        lines = f.readlines()
    
    count = 0 
    # print(lines[7:15])
    time_location = [] 
    for line in lines: 
        if line[2] == 'T':
            # print(count)
            time_location.append(count)
        count += 1 
        
    if len(time_location) == 1: 
        time_location.append(len(lines))
    
    # print(time_location)
    
    lines = lines[time_location[0]+7:time_location[1]-1]
    # print(lines[0:15])
    # print(lines[-1])
    
    atom_bond = {}
    C_ids = []
    for line in lines: 
        line = line.split(' ')
        line = list(filter(lambda item: item.strip(), line))
        line = [float(value) for value in line]
        atom_type = int(line[1]) 
        atom_id = int(line[0])
        # print(atom_type)
        # print(line)
        if atom_type == 2: 
            C_ids.append(atom_id)
        if atom_type == 3: 
            atom_bond_end_index = int(line[2]+3)
            atom_bond[atom_id] = line[3:atom_bond_end_index]
            
            
    Si_ids = list(atom_bond.keys())
    counterSi = 0
    counterC = 0
    for bonds in atom_bond.values(): 
        
        for bond in bonds:
            if bond in Si_ids:
                counterSi += .5
            if bond in C_ids:
                counterC += .5
    # print(counter)
        
    
    # print(len(Si_id),len(atom_bond))
    return [counterSi,counterC]
    

def Si_Csingle(file_name):
    with open(file_name) as f: 
        lines = f.readlines()
    
    count = 0 
    # print(lines[7:15])
    time_location = [] 
    for line in lines: 
        if line[2] == 'T':
            # print(count)
            time_location.append(count)
        count += 1 
        
    if len(time_location) == 1: 
        time_location.append(len(lines))
    
    lines = lines[time_location[0]+7:time_location[1]-1]
    # print(lines[0:15])
    # print(lines[-1])
    
    atom_bond = {}
    C_ids = []
    
    for line in lines: 
        line = line.split(' ')
        line = list(filter(lambda item: item.strip(), line))
        line = [float(value) for value in line]
        atom_type = int(line[1]) 
        atom_id = int(line[0])
        # print(atom_type)
        # print(line)
        if atom_type == 2: 
            C_ids.append(atom_id)
        if atom_type == 3: 
            atom_bond_end_index = int(line[2]+3)
            atom_bond[atom_id] = line[3:atom_bond_end_index]
            
            
    
    counter = 0
    for bonds in atom_bond.values(): 
        
        for bond in bonds:
            if bond in C_ids:
                counter += .5
    # print(counter)
        
    
    # print(len(Si_id),len(atom_bond))
    return counter
        
        
    
# Si_Sisingle('gas_creation.1.reaxc')
# Si_Sigroup()
Si_Simultidir()