# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 10:20:26 2021

@author: Harrison
"""

#specific gas deletions for certain atoms 
#MW_threshold does nothing 

def MW_gas_delete(lmpptr,bond_file,MW_threshold): 
    
    from lammps import lammps 
    lmp = lammps(ptr=lmpptr)
    
    with open('type_mass.csv','r') as f: 
        mass_types = f.readlines() 
    
    masses = {}
    for mass in mass_types: 
        mass = mass.split(',')
        atom_type = int(mass[0])
        atom_mass = float(mass[1])
        
        masses[atom_type] = atom_mass
        
    # print(masses)
        
    
    
    with open(bond_file) as f: 
        lines = f.readlines()

        
    split_locations = [] 
    x = 0 
    
    atom_info = {}
    
    for line in lines :
        if line[2] == 'T':
            split_locations.append(x)
        x += 1
        
    split_locations.append(len(lines))
        
    #print(split_locations)
    #print(len(lines))
    
    #atom id and mol id  {mol id, atom id}
    # mol_info = {}
    
    #atom id and id of atoms bonded to it 
    atom_bond = {}
    for line in lines[split_locations[-2] + 7:-1]:
        line = line.split(' ')
        line = list(filter(lambda item: item.strip(), line))
        line = [float(value) for value in line]
        atom_id = int(line[0]) 
        # mol_id_index = int(line[2]) + 3
        # mol_id = line[mol_id_index]
        # if mol_id not in mol_info: 
        #     mol_info[mol_id] = []
        #     mol_info[mol_id].append(atom_id)
        # else:
        #     mol_info[mol_id].append(atom_id)
        
        atom_info[atom_id] = []
        atom_info[atom_id].append(line[1])
        atom_mass = masses[atom_info[atom_id][0]]
        atom_info[atom_id].append(atom_mass)
        
        
        atom_bond_end_index = int(line[2]+3)
        atom_bond[atom_id] = (line[3:atom_bond_end_index])
            
        
    #molecule id 
            
    
    
    # print(atom_bond[1536])
    # print(timesteps_refined[0][0:10])
    # return(timesteps_refined)
    
    # print(split_locations)
    # print(mol_info)
    
    total_atom = len(atom_bond)
    atom_check_list = [*range(1,total_atom+1)]
    # print(atom_check_list)
    
    mol_num = 1 
    mol_info = {}
    
    
    #algorithm to 
    while len(atom_check_list) != 0:
        
        mol_info[mol_num] = []
        temp_mol_check = [atom_check_list[0]]
        atom_check_list.remove(temp_mol_check[0])
    
        
        # print(len(mol_info[mol_num]))
        
        
        while len(temp_mol_check) != 0: 
            
    
            atom_bonding = atom_bond[temp_mol_check[0]]
            # print(atom_bonding)
            
            for bonds in atom_bonding:
                if bonds not in temp_mol_check and bonds in atom_check_list:
                    temp_mol_check.append(bonds)
                    atom_check_list.remove(bonds)
                    
            # print(len(temp_mol_check))
            # print(atom_check_list)
            # print(temp_mol_check)
            
            
            x = temp_mol_check[0]
            mol_info[mol_num].append(x)
            temp_mol_check.remove(temp_mol_check[0])
            
            
            # if x == y: 
            #     print(temp_mol_check)
            #     temp_mol_check =[]
                    
            
        
        mol_num = mol_num + 1
        
    # print(mol_info[6])
    
    #mol_mass = {id:[mass,Si(0 if no)}
    mol_mass = {}
    for values in [*range(1,mol_num)]: 
        mol_mass[values] = [0,0]
    # print(mol_mass)
    
    # print(mol_mass)
    
    for values in range(1,mol_num):
        sample_mass = 0 
        Si_counter = 0 
        for ids in mol_info[values]:
            sample_mass = sample_mass + atom_info[ids][1]
            if atom_info[ids][1] >= 25: 
                Si_counter += 1 
                
        # print(sample_mass)
        # print(Si_counter)
        # print(mol_mass[values][0])
        mol_mass[values][0] = sample_mass
        mol_mass[values][1] = Si_counter
        
    
    
    deletion_list = [16.04246, 30.06904, 28.05316, 26.03728, 2.01588, 28.0101, 44.0095, 18.01528, 31.9988, 78.11184]
    mol_del = []
            
    for mass_keys in mol_mass.keys(): 
        if round(mol_mass[mass_keys][0],5) in deletion_list: 
            # print(mol_mass[mass_keys])
            mol_del.append(mass_keys)
    
    atom_del = []
    for values in mol_del: 
        # print(mol_info[values])
        atom_del = atom_del + mol_info[values]
        
    # print(atom_del)
    
    if len(atom_del) > 0: 
        atom_del_string =''
        
        for atom_id in atom_del: 
            atom_del_string = atom_del_string + str(int(atom_id)) + ' '
            
        group_str = "group gas id " + atom_del_string

        lmp.command(group_str)
        lmp.command("delete_atoms group gas")
        print("gases deleted")
        
    
MW_gas_delete(7,'gas_creation.15.reaxc',40)
