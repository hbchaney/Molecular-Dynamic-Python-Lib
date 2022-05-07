# -*- coding: utf-8 -*-
"""
Created on Sat Jan 29 13:29:16 2022

@author: hbcha
"""

import os
import json 
import matplotlib.pyplot as plt
from copy import deepcopy


def Bondmultdirect(parent_dir,atom_types,data_out,name):
    '''
    

    Parameters
    ----------
    parent_dir : string
        location of all the directories holding bond info 
    atom_types : list of integers
        [Si,O,C,H]
    data_out : string
        sends the data to a direct in this location
    name : 
        name of output file 

    Returns
    -------
    a json file of all bond info 

    '''
    
    Si = atom_types[0]
    O = atom_types[1]
    C = atom_types[2]
    H = atom_types[3]
    
    os.chdir(parent_dir)
    
    dir_present = sorted(os.listdir())
    # print(dir_present)
    cwd = os.getcwd()
    
    data = []
    count = 0 
    t = False
    
    
    
    for directs in dir_present: 
        os.chdir(directs)
        
        data_new = Bondsingledirect(Si,O,C,H,t)
        
        
        for points in data_new: 
            points['other']['timestep'] += count
            
        data += data_new
        
        count = data_new[-1]['other']['timestep']
        
        
        os.chdir(cwd)
        
        t = True
        
    with open(data_out+'/'+name,'w') as f:
        json.dump(data,f,indent=3)
        
    print('done file created')

    

def Bondsingledirect(Si,O,C,H,t=True):
    '''
    

    Parameters
    ----------
    Si : int
        int corresponding with Si 
    O : int
        int corresponding with O
    C : int
        int corresponding with C 
    H : int
        int corresponding with H 

    Returns
    -------
    a list of dictionaries with a corresponding list [{Si:[Si,O,C,H],O:[Si,...]},{Si:[Si,O,...]}

                                                     
    '''
    
    type_coord = {Si:'Si',O:'O',C:'C',H:'H'}
    bond_coord = {Si:0,O:1,C:2,H:3}
    
    files = os.listdir()
    
    files = sorted(files, key=lambda x: float(x[13:-6]))
    
    # print(files)
    
    
    data = []
    
    if t != True :
        data += Bond_singlefile_first(bond_coord,type_coord,files[0])
    for file in files: 
        data_new = Bond_singlefile(bond_coord,type_coord,file)
        
        # for values in data_new : 
        #     print(values['other']['timestep'])
        data += data_new
        
        
    # print(data)
    # for dics in data:
    #     print(dics['Si'][2])
                
            
    # print(data)
    return data


def Bond_singlefile_first(bond_coord,type_coord,file_name) : 
    
    bond = {'Si':[0,0,0,0],'O':[0,0,0,0],'C':[0,0,0,0],'H':[0,0,0,0],'other':{'num':[0,0,0,0],'timestep':0}}
    
    with open(file_name,'r') as f:
        lines = f.readlines()
        
    
    data = []
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
    # print(lines[0:30])
    
    atom_id_type = {}
    refined_list = []
    for line in lines:
        line = line.split(' ')
        line = list(filter(lambda item: item.strip(), line))
        line = [float(value) for value in line]
        refined_list.append(line)
        # print(line)
        
        atomtype = bond_coord[line[1]]
        
        bond['other']['num'][atomtype] += 1 
        
        atom_id_type[line[0]] = line[1]
        # print(atom_id_type)
        
    for line in refined_list: 
        
        atom_id = line[0]
        atom_type = type_coord[atom_id_type[atom_id]]
        # print(atom_type)
        
        #bonds 
        num_bonds = int(line[2])
        bond_ids = line[3:(num_bonds+3)]
        
        for ids in bond_ids: 
            
            bond_type = bond_coord[atom_id_type[ids]]
            bond[atom_type][bond_type] += 1 
            
    # print(bond)
    data.append(bond)
    return data 
    
    
def Bond_singlefile(bond_coord,type_coord,file_name): 
    
    
    with open(file_name,'r') as f:
        lines = f.readlines()
        
    
    data = []
    count = 0 
    # print(lines[7:15])
    time_location = [] 
    for line in lines: 
        if line[2] == 'T':
            # print(count)
            time_location.append(count)
        count += 1 
        
    end = len(lines) #because not inclusive 
    time_location.append(end)
    
    # print(time_location)
    
    if len(time_location) == 2:
        return []
    
        
    
    count = 1 
    while count < len(time_location)-1 : 
        
        bond = {'Si':[0,0,0,0],'O':[0,0,0,0],'C':[0,0,0,0],'H':[0,0,0,0],'other':{'num':[0,0,0,0],'timestep':0}}
        
        lines_temp = lines[time_location[count]+7:time_location[count+1]-1]
        # print(lines[0:30])
        
        # print(lines[time_location[count]])
        
        ts = int(lines[time_location[count]].split()[2])
        # print(ts)
        bond['other']['timestep'] = ts
        
        atom_id_type = {}
        refined_list = []
        for line in lines_temp:
            line = line.split(' ')
            line = list(filter(lambda item: item.strip(), line))
            line = [float(value) for value in line]
            refined_list.append(line)
            # print(line)
            
            atomtype = bond_coord[line[1]]
            
            bond['other']['num'][atomtype] += 1 
            
            atom_id_type[line[0]] = line[1]
            # print(atom_id_type)
            
        for line in refined_list: 
            
            atom_id = line[0]
            atom_type = type_coord[atom_id_type[atom_id]]
            # print(atom_type)
            
            #bonds 
            num_bonds = int(line[2])
            bond_ids = line[3:(num_bonds+3)]
            
            for ids in bond_ids: 
                
                bond_type = bond_coord[atom_id_type[ids]]
                bond[atom_type][bond_type] += 1 
                
        # print(bond)
        data.append(bond)
        
        count = count + 1
        
        
    return data
        
    
    
    
    
    
    
def json_plot(json_file):
    
    with open(json_file,'r') as f:
        contents = json.load(f)
        
    x = []
        
    Si_C = []
    Si_H = []
    Si_Si = []
    Si_O = []
    Si_all = []
    C_Si = []
    num_car = []
    num_sil = []
    num_hydro = []
    num_oxy = []
    C_H = []
    C_C = []
    
    wt_car = []
    wt_sil = []
    wt_hydro = []
    wt_oxy = []
    
    total = []
    tatoms = []
    
    starting_weight = contents[0]['other']['num'][2]*12+contents[0]['other']['num'][0]*28+contents[0]['other']['num'][1]*16+contents[0]['other']['num'][3]
    start = contents[0]
    start_C_H = start['C'][3]/start['other']['num'][2]
    
    
    for dicts in contents:
        x.append(dicts['other']['timestep']*.2)
        Si_C.append(dicts['Si'][2]/dicts['other']['num'][0]) 
        Si_H.append(dicts['Si'][3])
        Si_Si.append(dicts['Si'][0]/2/dicts['other']['num'][0])
        Si_O.append(dicts['Si'][1])
        
        C_H.append(dicts['C'][3]/dicts['other']['num'][2]/start_C_H)
        C_C.append(dicts['C'][2]/2/dicts['other']['num'][2])
        
        num_c = dicts['other']['num'][2]
        num_car.append(num_c)
        
        num_sil.append(dicts['other']['num'][0])
        num_oxy.append(dicts['other']['num'][1])
        num_hydro.append(dicts['other']['num'][3])
        
        t_atm = dicts['other']['num'][0] + dicts['other']['num'][1] + dicts['other']['num'][3] + dicts['other']['num'][2]
        
        t1= (num_c*12+dicts['other']['num'][0]*28+dicts['other']['num'][1]*16+dicts['other']['num'][3])/starting_weight
        t= (num_c*12+dicts['other']['num'][0]*28+dicts['other']['num'][1]*16+dicts['other']['num'][3])

        total.append(t1)
        tatoms.append(t_atm)
        
        wt_car.append(num_c*12/t*100)
        wt_sil.append(dicts['other']['num'][0]*28/t)
        wt_hydro.append(dicts['other']['num'][3]/t)
        wt_oxy.append(dicts['other']['num'][1]*16/t)
        
        C_Si.append(dicts['C'][0]/num_c)
        # C_Si.append(dicts['C'][0])
        
    
    # print(len(Si_C))
        
    #t8   
    # x = list(range(0,len(Si_C)*5000,5000))
    
    # t9 
    # x = list(range(0,40*5000,5000)) + list(range(40*5000,((len(Si_C)-40)*5000*2+40*5000),5000*2))
    
    #t10 
    # x = list(range(0,80*5000,5000)) + list(range(80*5000,((len(Si_C)-80)*5000*2+80*5000),5000*2))

    #t11 
    # x = list(range(0,120*5000,5000)) + list(range(120*5000,((len(Si_C)-120)*5000*2+120*5000),5000*2))

    
    # print(x)
    # plt.plot(x,Si_C,label = 'Si-C')
    # plt.plot(x,C_H,label = 'C-H')
    # plt.plot(x,C_C,label='C-C')
    # plt.plot(x,Si_H,label = 'Si-H')
    # plt.plot(x,Si_Si,label = 'Si-Si')
    # plt.plot(x,Si_O,label = 'Si-O')
    # plt.plot(x,Si_all,label = '% bonded')
    # plt.plot(x,C_Si,label = 'C-Si bonds ')
    # plt.plot(x,num_car,label = 'number of carbons')
    plt.plot(x,num_oxy,label='oxygen')
    # plt.plot(x,num_sil,label='silicon')
    # plt.plot(x,num_hydro,label='Hydrogen')
    # plt.plot(x,wt_sil,label = 'Si')
    # plt.plot(x,wt_hydro,label = 'H')
    # plt.plot(x,wt_oxy,label = 'O')
    # plt.plot(x,wt_car,label = 'C%')
    # plt.plot(x,total,label = 'wt%')
    # plt.plot(x,tatoms,label = 'total atoms')
    
    # plt.title("change in C-H bonds as a fraction of starting C-H bonds")
    # plt.xlabel('time in femtoseconds')
    # plt.ylabel('number of C-H bonds as fraction of start')
    
    # plt.ylim(ymin=0)
    plt.legend()
    
    # plt.show()
    
def json_bar(json_file):
    
    with open(json_file,'r') as f:
        contents = json.load(f)
        
    intial = contents[-1]['other']['num']
    
    Si, O, C, H = intial 
    
    total = Si*28 + O*16 + C*12.01 + H 
    
    percents = [Si*28/total,O*16/total,C*12/total,H/total]
    names = ['Si','O','C','H']
    
    fig, ax = plt.subplots()
    ax.bar(names,percents)
    
    for i, v in enumerate(percents):
        print(i,v)
        plt.text(x = i-.1 , y = v , s = f'{round(v,3)}',fontdict=dict(fontsize=10))
        
    
    print(Si,O,C,H)
        
        



# Bondmultdirect(r'G:\My Drive\Research\reaxff\PSO\PSO_t12_gas',[1,2,3,4],r'G:\My Drive\Research\reaxff\PSO\gas_results','PSO_t12_gas.json')
# Bondmultdirect('G:/My Drive/Research/reaxff/mc/t7_gas_data', [3,4,2,1], 'G:/My Drive/Research/reaxff/mc/gas_results','mc-t7-gas.json')
# Bondmultdirect('G:/My Drive/Research/reaxff/mc/t5_gas_data',[3,4,2,1],'G:/My Drive/Research/reaxff/mc/gas_results','mc-t5-gas.json')
# Bondmultdirect('G:/My Drive/Research/reaxff/mc/PSO_t8_gas', [3,4,2,1], 'G:/My Drive/Research/reaxff/mc/gas_results', 'mc-t8-gas.json')
# Bondmultdirect(r'G:\My Drive\Research\reaxff\PSO\PSO_t11_gas', [3,4,2,1], r'G:\My Drive\Research\reaxff\PSO\gas_results', 'mc-t11-gas.json')
# Bondmultdirect(r'G:\My Drive\Research\reaxff\PHMS\PHMS_gas', [4,3,2,1], r'G:\My Drive\Research\reaxff\PHMS\PHMS_gas_results', 'PHMS_t1_1.json')
# Bondmultdirect(r'G:\My Drive\Research\reaxff\VMHM\VMHM_t2_gas', [1,2,3,4], r'G:\My Drive\Research\reaxff\VMHM\VMHM_gas_results', 'VMHM_t2.json')
# json_plot(r'G:\My Drive\Research\reaxff\PSO\gas_results\mc-t11-gas.json')
# Bondmultdirect(r'G:\My Drive\Research\reaxff\PDVS\PDVS_t2_gas',[1,2,3,4],r'G:\My Drive\Research\reaxff\PDVS\PDVS_gas_data','PDVS_gas_data.json')
# json_plot(r'G:\My Drive\Research\reaxff\PDVS\PDVS_gas_data\PDVS_gas_data.json')

# Bondmultdirect(r'G:\My Drive\Research\reaxff\PDMS\PDMS_t2_gas', [1,2,3,4], r'G:\My Drive\Research\reaxff\PDMS\PDMS_gas_data', 'PDMS_t2_gas.json')
# json_plot(r'G:\My Drive\Research\reaxff\PHMS\PHMS_gas_results\PHMS_t1_1.json')
# json_plot(r'G:\My Drive\Research\reaxff\VMHM\VMHM_gas_results\VMHM_t2.json')
# json_bar(r'G:\My Drive\Research\reaxff\PDMS\PDMS_gas_data\PDMS_t2_gas.json')

json_plot(r'G:\My Drive\Research\reaxff\PSO\gas_results\PSO_t12_gas.json')
