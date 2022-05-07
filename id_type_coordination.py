# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 08:23:38 2021

@author: Harrison
"""

def infile(file): 
    '''
    reads the input file and creates an easy to parse csv with the 
    id,type,mass data in each line 
          
    '''
    
    
    
    with open(file,'r') as f: 
        lines = f.readlines() 
        
    #printing out seperate catagories 
    
    #Masses 
    # print(lines[9:13]) #the masses 
    
    masses = []
    
    for values in lines[9:13]:
        mass_info = values.split(' ')
        # print(mass_info)
        
        atom_type = int(mass_info[0])
        atom_mass = float(mass_info[1])
        
        masses.append([atom_type,atom_mass])
        
    print(masses)
    
    #finding where the info is stored for the atom ids 
    # print(lines[16:20])
    # print(atom_info)
    
    with open('type_mass.csv','w') as f:
        for values in masses: 
            f.write(f'{values[0]},{values[1]}\n')
        print('atom_mass info csv created')
    
    

infile('rand_2x2_t4_500000.in')