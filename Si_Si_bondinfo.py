# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 14:01:31 2021

@author: Harrison
"""

import os 
import matplotlib.pyplot as plt
# print(os.listdir())
os.chdir('stable_pryo_t4_bondparse')
# print(os.listdir())

def Si_Sigroup(): 
     files = os.listdir()
     # print(sorted(files))
     # print(files[0][13:-6])
     files = sorted(files, key=lambda x: float(x[13:-6]))
     # print(files)
     ts = range(0,len(files))
     
     Si_Si = []
     
     for file in files: 
         bonds = Si_Sisingle(file)
         Si_Si.append(bonds)
         
     # print(Si_Si)
     plt.plot(ts,Si_Si)
     plt.xlabel('timesteps')
     plt.ylabel('Si-Si bonds')
     plt.show()
     
     
     
     


    
    
# 'gas_creation.1.reaxc'   
def Si_Sisingle(file_name): 
    
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
    
    lines = lines[time_location[0]+7:time_location[1]-1]
    # print(lines[0:15])
    # print(lines[-1])
    
    atom_bond = []
    Si_id = []
    for line in lines: 
        line = line.split(' ')
        line = list(filter(lambda item: item.strip(), line))
        line = [float(value) for value in line]
        atom_type = int(line[1]) 
        atom_id = int(line[0])
        # print(atom_type)
        # print(line)
        if atom_type == 3: 
            Si_id.append(atom_id)
            atom_bond_end_index = int(line[2]+3)
            atom_bond.append(line[3:atom_bond_end_index])
            
    
    
    # print(len(Si_id),len(atom_bond))
    # print(atom_bond[0:15])
    total_bonds = 0
    for bonds in atom_bond: 
        total_bonds += len(bonds) 
        
    # print(total_bonds/len(atom_bond))
        
    Si_Si_count = 0 
    for bonds in atom_bond: 
        # print(bonds)
        if type(bonds) == list:
            for bond in bonds: 
                # print(bond)
                if bond in Si_id: 
                    Si_Si_count += 1 
        else: 
            if bonds in Si_id: 
                Si_Si_count += 1
    
    Si_Si_count = Si_Si_count/2
    # print(Si_Si_count)
    return Si_Si_count
        
        
    
# Si_Sisingle('gas_creation.25.reaxc')
Si_Sigroup()