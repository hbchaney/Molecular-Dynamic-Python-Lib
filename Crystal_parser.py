# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 13:51:27 2022

@author: hbcha
"""

import numpy as np
from copy import deepcopy 
import datetime

class bond_frame:
    '''
    gas_frame data structure 
                    [id]     [0]    [1]         [2][0]   
    bond_data = {atm_id:[atm_type,num_bonds,[atm_ids_bonded_to]],atm_id2:...}
    timestep = int value (used as a check with loc_frame) (also might be used 
                                                           in bond info tracking)
    '''
    
    def __init__(self,bonding=[],ts=0):
        self.bond_data = bonding 
        self.num_atoms = len(bonding) 
        self.ts = ts
        
    def duplicate(self): 
        
        new = bond_frame()
        new.bond_data = deepcopy(self.bond_data)
        new.num_atoms = self.num_atoms 
        new.ts= self.ts 
        
        return new
        
        
    def file_read(self,file_name,coord=[1,2,3,4]): 
        '''

        Parameters
        ----------
        cls : bond_frame or subclass with bond_data attr
        file_name : str
            name of file where the bonding data is stored

        Returns
        -------
        None.

        '''
        type_coord = {coord[0]:1,coord[1]:2,coord[2]:3,coord[3]:4}
        
        
        with open(file_name, 'r') as f: 
            data = f.read() 
           
        #parsing the file for last timestep 
        lines = data.split('\n') 
        i = 0 
        ts_loc = [] #list of ts indicies 
        for values in lines:
            if values[0:3] == '# T' :
                ts_loc.append(i)
                 
            i += 1
            
        # setting self.ts = to the last ts value 
        ts_str = lines[ts_loc[-1]]
        ts = int(ts_str.split()[-1]) #ts set 
        
        #start and end of data in question 
        start = ts_loc[-1] + 7
        end = len(lines)-2
        
        #relevent info stored in atom_info
        atom_info = lines[start:end] 
        
        #init bond_data
        bond_data = {}
        
        for info in atom_info:
            #getting rid of empty spaces 
            info = info.split(' ')
            info = ' '.join(info).split()
            
            atm_id = int(info[0]) 
            #converting the atm_type based on coord to std coord (1,2,3,4)
            atm_type = type_coord[int(info[1])]
            bond_num = int(info[2])
            bonded_atms = info[3:3+bond_num]
            
            bonded_atms = [int(values) for values in bonded_atms] 
            
            # print([atm_type,bond_num,bonded_atms])
            bond_data[atm_id] = [atm_type,bond_num,bonded_atms]
        
        self.bond_data = bond_data 
        self.ts = ts 
        self.num_atoms = len(bond_data) 
            
            
        
    def bond_exclude(self,atm_type1,atm_type2,reverse=0) : 
        '''
        Parameters
        ----------
        atm_type1 : int
            first atom type 
        atm_type2 : int
            second atom type 
        
        uses std coordination numbers 
        goes through every atom in bond_data and checks the if the 
        atm_type matches atm_type 1 or 2 it checks to see if it is bonded to 2 
        it then saves the atm_id to be deleted from loc_data 
        
        Returns
        -------
        list of atm_ids to be deleted 
        '''
        
        
        bond_list = [] #list of ids 
        bond_rlist = []
        count = 0
        for k,v in self.bond_data.items() :
            #checks to see if atm_type for k is part of atm_type1 or atm_type2 
            bo = 0
            if v[0] == atm_type1 : 
                for atm_bonds in v[2]: 
                    if atm_bonds not in self.bond_data: 
                        count += 1
                    elif self.bond_data[atm_bonds][0] == atm_type2 :
                        bo = 1
            if v[0] == atm_type2 : 
                for atm_bonds in v[2]: 
                    if atm_bonds not in self.bond_data: 
                        count += 1
                    elif self.bond_data[atm_bonds][0] == atm_type1 :
                        bo = 1
            
            if bo == 1 : 
                bond_list.append(k)
                
            else : 
                bond_rlist.append(k)
                
        print(count)
                
        if reverse == 0:      
            return bond_list  
        else : 
            return bond_rlist 
    
    
    def _bond_del(self,del_list) :
        '''

        Parameters
        ----------
        del_list : list of ints
            deletes all the atoms in atm list and gets rid of the bonds that 

        Returns
        -------
        None.

        '''

        for values in del_list: 
            del self.bond_data[values] 
        
        self.num_atoms = len(self.bond_data) 

        
        for k,v in self.bond_data.items(): 
            for bnd_numbers in self.bond_data[k][2] : 
                if bnd_numbers not in self.bond_data.keys(): 
                    self.bond_data[k][2].remove(bnd_numbers)

        
            
            
    
            
class loc_frame: 
    '''
    holds the locational data with the atom ids 
    loc_data = {atom_id:[atm_type,(x,y,z)]}
    timestep = int value (used as a check with bond data)
    '''
    
    def __init__(self,locations=[],ts=0) :
        self.loc_data = locations
        self.num_atoms = len(locations) 
        self.ts = ts
        self.box = [] 
        
    def duplicate(self) : 
        new = loc_frame() 
        new.loc_data = deepcopy(self.loc_data) 
        new.box = deepcopy(self.box) 
        new.ts = self.ts 
        new.num_atoms = self.num_atoms
        
        return new
    def to_file(self,file_name): 
        '''
        Parameters
        ----------
        file_name : str
            name and full path to file
            
        Takes in the data from loc_dat and creates a file from it.  
        the atom ids are lost because they are not expected to have any meaning 

        Returns
        -------
        None.

        '''
        
        
        num_atoms = len(self.loc_data)
        
        masses = [28.0855,15.9994,12.0107,1.00794]
        
        atms = []
        for values in self.loc_data.values(): 
            atms.append([values[0],values[1][0],values[1][1],values[1][2]])
        atms = np.array(atms)
            
        # print(atms.T)
        Master = []
        
        Master.append(f'LAMMPS data file written by Harrison made on {datetime.date.today()}')
        Master.append(f'{num_atoms} atoms')
        Master.append('4 atom types')
        
        xhi,xlo = np.amax(atms.T[1]),np.amin(atms.T[1])
        yhi,ylo = np.amax(atms.T[2]),np.amin(atms.T[2])
        zhi,zlo = np.amax(atms.T[3]),np.amin(atms.T[3])
        # print(xlo,xhi,ylo,yhi,zlo,zhi)
            
        Master.append(f'{xlo} {xhi} xlo xhi')
        Master.append(f'{ylo} {yhi} ylo yhi')
        Master.append(f'{zlo} {zhi} zlo zhi')
        
        Master.append('')
        Master.append('Masses')
        Master.append('')
        
        Master.append(f'1 {masses[0]} #Si')
        Master.append(f'2 {masses[1]} #O')
        Master.append(f'3 {masses[2]} #C')
        Master.append(f'4 {masses[3]} #H')
        
        Master.append('')
        Master.append('Atoms  # charge')
        Master.append('')
        
        i = 1
        for values in atms: 
            Master.append(f'{i} {int(values[0])} 0.0 {values[1]} {values[2]} {values[3]}')
            i += 1
        Master_line =[]
        
        for values in Master:
            values = values + '\n'
            Master_line.append(values)
            
        Master = Master_line
        
            
        with open(file_name,'w') as f:
            f.writelines(Master)
            
        print('done')
        
    def file_read(self,file_name,coord=[1,2,3,4]): 
        
        '''
        Parameters
        ----------
        cls : loc_frame or subclass with loc_data variable
        file_name : str
            name of file where the location data is stored and filename 

        Returns
        -------
        None.

        '''
        type_coord = {coord[0]:1,coord[1]:2,coord[2]:3,coord[3]:4}
        
        with open(file_name, 'r') as f: 
            data = f.read()
        
        lines = data.split('\n') 
        #finding last time step 
        
        ts_loc = []
        i=0 
        for values in lines: 
            if values[0:7] == 'ITEM: T': 
                ts_loc.append(i+1)
            i += 1
            
        #box creation 
        box_s = ts_loc[-1] + 4
        box_e = ts_loc[-1] + 7 
        box = lines[box_s:box_e]
        temp_list = []
        for things in box :
            temp = things.split(' ')
            temp = [float(val) for val in temp]
            temp_list.append(temp)
            
        box = np.array([temp_list[0],temp_list[1],temp_list[2]])
        ran = box.T[1] - box.T[0]
        
        self.box = box 
        
        #timestep     
        ts = int(lines[ts_loc[-1]])
        
        #making the start and end 
        start = ts_loc[-1] + 8
        end = len(lines)-1
        
        #relevent data 
        loc_data = {}
        atom_data = lines[start:end] 
        
        for values in atom_data:
            values = values.split() 
            atm_id = int(values[0])
            atm_type = type_coord[int(values[1])]
            xyz = np.array([float(values[2]),float(values[3]),float(values[4])])
            xyz = xyz * ran 
            xyz = xyz + box.T[0]
            
            loc_data[atm_id] = [atm_type,(xyz[0],xyz[1],xyz[2])]
            
        # print(loc_data)
        self.loc_data = loc_data 
        self.num_atoms = len(loc_data) 
        self.ts = ts 
        
    def _loc_del (self,id_list) : 
        '''

        Parameters
        ----------
        id_list : list of int
            takes in list of int and del the 

        Returns
        -------
        None.

        '''
        for values in id_list : 
            del self.loc_data[values] 
            
        self.num_atoms = len(self.loc_data)
        
    def __str__(self): 
        atm_types = {1:0,2:0,3:0,4:0} 
        for values in self.loc_data.values() : 
            atm_types[values[0]] += 1 
        Si = f'{atm_types[1]}'
        O = f'{atm_types[2]}'
        C = f'{atm_types[3]}'
        H = f'{atm_types[4]}'
        
        tot = f'# of atoms:\nSi : {Si}\nO : {O}\nC : {C}\nH : {H}\ntot = {int(Si) + int(O) + int(C) + int(H)}'
        return tot
            
            
            
        
    
    
class combined_frame(bond_frame,loc_frame): 
    
    '''
    sub class of both bond frame and loc_frame. 
    used to identify and delete atoms to look at certain regions in isolation 
    '''
    
    def __init__(self,bonding=[],locations=[]) : 
        self.loc_data = locations 
        self.bond_data = bonding 
        self.num_atoms = len(locations) 
        self.box = []
        self.ts = 0 
        
    def duplicate(self): 
        new = combined_frame()
        new.loc_data = deepcopy(self.loc_data) 
        new.box = deepcopy(self.box) 
        new.bond_data = deepcopy(self.bond_data)
        new.num_atoms = self.num_atoms 
        new.ts= self.ts 
        
        return new
        
        

    def file_read(self,gas_file,loc_file,coord=[1,2,3,4]) : 
        loc_frame.file_read(self,loc_file,coord=coord)
        loc_ts = self.ts 
        loc_num = self.num_atoms
        bond_frame.file_read(self,gas_file,coord=coord) 
        
        if loc_ts != self.ts: 
            raise ValueError('timesteps not consistent')
        if loc_num != self.num_atoms: 
            raise ValueError('number of atoms is not consistent')
            
            
    def comb_del(self,del_list): 
        self._bond_del(del_list)
        self._loc_del(del_list) 
        
    
        
        
            
    
        
    
        
        

    
        
    
    
    