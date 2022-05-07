# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 13:25:56 2022

@author: Harrison
"""

'''
the intention of this file is to be called upon as a module 

the format of atom info stored will no longer have the atom id 

it will be as follows [(x,y,z,'c'),(x,y,z,'Si'),...]
'''

import numpy as np 
import math 
from random import randint
import datetime
import re 


#####-------------------------------------------------------------------##### 
#####-----------------------------CONSTANTS-----------------------------#####
#####-------------------------------------------------------------------#####
    
class Bond_Length: 

    #all values in angstroms 
    C_CPh = 1.39 #carbon in ring 
    C_C = 1.535 
    C_Cd = 1.339 #carbon double bond 
    C_H = 1.094 
    C_O = 1.43
    #Si_O = 1.483 increased the length of the bond 
    Si_O = 1.583
    Si_C = 1.86 
    Si_H = 1.58 
    
#####-------------------------------------------------------------------##### 
#####-----------------------------TOOLS---------------------------------#####
#####-------------------------------------------------------------------#####

def Rotation_func(coords_in,psi,theta,phi): 

    '''
    takes in a tuple containing the coordinates in tuple form, and the atom type
    
    all rotations are clockwise?
    all rotation inputs in degrees 
    
    psi = rotation about the x axis 
    theta = rotation about the y axis 
    phi = rotation about the z axis 
    
    '''
    
    psi,theta,phi = np.radians(np.array([psi,theta,phi]))
    
    #creating the matricies of rotation 
    psi_matrix = np.array([[1, 0, 0], 
                           [0, np.cos(psi), np.sin(psi)],
                           [0, -1*np.sin(psi), np.cos(psi)]])
    
    theta_matrix = np.array([[np.cos(theta), 0, np.sin(theta)],
                              [0, 1, 0],
                              [-1*np.sin(theta), 0, np.cos(theta)]])
    
    phi_matrix = np.array([[np.cos(phi), np.sin(phi), 0],
                            [-1 *np.sin(phi), np.cos(phi), 0], 
                            [0, 0, 1]])
                            
    #combining the matricies to make a general rotation matrix 
    combined_matrix = phi_matrix @ theta_matrix @ psi_matrix
    # print(combined_matrix)
    
    new_list = []
    
    for coord in coords_in : 
        coords  = np.array(coord[0:3])
        # print(coords)
        El = coord[3]
        
        new_c = (coords @ combined_matrix)
        # print(new_c)
        
        new_list.append((new_c[0],new_c[1],new_c[2],El))
        
    return new_list
        

def Translation_func(coord_in,trans): 
    
    '''
    takes in a list of coordinates and and specific translation and outputs 
    the a list translated coordiantes 
    
    can also output individual coordinates if needed 
    
    trans input: 
        [x,y,z] or (x,y,z)
    '''
    
    trans = np.array(trans)
    new_list = []
    
    for coord in coord_in:
        coords = np.array(coord[0:3])
        El = coord[3]
        
        new_c = coords + trans
        
        new_list.append((new_c[0],new_c[1],new_c[2],El))
        
    return new_list
    

def Mirror_func(coord_in,plane): 
    '''
    takes in a list of coordinates and reflects them over a specific plane 
    
    input plane:
        'x' 'y' 'z'
    '''
    
    x_m = np.array([-1,1,1])
    y_m = np.array([1,-1,1])
    z_m = np.array([1,1,-1])
    
    if plane == 'x': 
        plane = x_m
    elif plane == 'y': 
        plane = y_m 
    elif plane == 'z':
        plane = z_m
        
    new_list = []
    
    for coord in coord_in: 
        coords = np.array(coord[0:3])
        El = coord[3]
        
        new_c = coords * plane
        
        new_list.append((new_c[0],new_c[1],new_c[2],El))
        
    return new_list
  

def file_maker(coordinates,out_loc,out_type = 'lammps'):
    '''
    

    Parameters
    ----------
    coordinates : list of tuples 
        
    out_type : 'lammps' or 'xyz' default 'lammps'
    
    out_loc : 'file_location and names'
    
    Returns
    -------
    prints if successful and the file type as well as the location
    makes a new file for later use 

    '''
    
    if out_type == 'lammps': 
        
        num_atoms = len(coordinates)
        
        masses = [28.0855,15.9994,12.0107,1.00794]
        
        Master = []
        
        Master.append(f'LAMMPS data file written by Harrison made on {datetime.date.today()}')
        Master.append(f'{num_atoms} atoms')
        Master.append('4 atom types')
        
        coord_str = []
        
        xlo,xhi = coordinates[0][0],coordinates[0][0]
        ylo,yhi = coordinates[0][1],coordinates[0][1]
        zlo,zhi = coordinates[0][2],coordinates[0][2]   
        
        type_coord = {'Si':1,'O':2,'C':3,'H':4}
        
        count = 1
        
        for coords in coordinates: 
            
            if coords[0] > xhi:
                xhi = coords[0]
            elif coords[0] <xlo: 
                xlo = coords[0]
                
            if coords[1] > yhi:
                yhi = coords[1]
            elif coords[1] < ylo:
                ylo = coords[1]
                
            if coords[2] > zhi:
                zhi = coords[2]
            elif coords[2] < zlo:
                zlo = coords[2]
                
                
            str_unit = f'{count} {type_coord[coords[3]]} 0.0 {round(coords[0],6)} {round(coords[1],6)} {round(coords[2],6)}'
            coord_str.append(str_unit)
            
            count += 1
            
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
        
        Master = Master + coord_str
        Master_line =[]
        
        for values in Master:
            values = values + '\n'
            Master_line.append(values)
            
        Master = Master_line
        
        # print(Master)
        
        with open(out_loc,'w') as f:
            f.writelines(Master)
            
        print('done')
        
    

#####-------------------------------------------------------------------##### 
#####-----------------------------SEGMENTS------------------------------#####
#####-------------------------------------------------------------------#####

def Benzene(): 
    
    '''
    Creates a benzene molecule coordinate points and outputs the values as a list 
    with the first value being the El symbol 
    
                       Hplace holder 
                       |
                       |
                      C1
                  ___/  \___
        H7----C2/           \C3-----H8
               |              |
               |              |
               |              |
        H9----C4             C5-----H10
                ----\   /----
                     C6
                     |
                     |
                     H11
    
    Starts at C3 and needs H3 to be present to work 
    any other H's can be commeneted out 
    will not output the specfic name of the atoms but rather just the El symbol 
    as the first value in the list 
    
    note: uses both the rotation function 
    
    '''
    
    #Initializing the position of C3 and H3 
    
    C3 = ((.5 * Bond_Length.C_CPh)/(math.tan(math.radians(30))), .5 * Bond_Length.C_CPh, 0, 'C')
    H8 = ((C3[0] + math.cos(math.radians(30))), C3[1] + math.sin(math.radians(30)), 0, 'H') 
    
    # C3, H8  = Rotation_func([C3,H8],0,0,0)
    C5, H10 = Rotation_func([C3,H8],0,0,60)
    C6, H11 = Rotation_func([C3,H8],0,0,120)
    C4, H9  = Rotation_func([C3,H8],0,0,180)
    C2, H7  = Rotation_func([C3,H8],0,0,240)
    C1      = Rotation_func([C3],0,0,300)[0]
    
    # H12 = Rotation(H8,0,0,300) #take out if needed 
    # H12[4] = 12
    
    return [C1,C2,C3,C4,C5,C6,H7,H8,H9,H10,H11] #take out H12
    

def Methal(): 
    
    '''
    creates a methal group by using an initial point and rotating around the main bonded axis 
    
    the carbon created exists at the point [0,0,0]
    
        
    '''
    C1 = [(0,0,0,'C')]
    
    H2 = [(math.cos(math.radians(19.5))*Bond_Length.C_H,
          -1 * math.sin(math.radians(19.5)) * Bond_Length.C_H,
          0,
          'H')]
    
    H3 = Rotation_func(H2,0,120,0)
    H4 = Rotation_func(H2,0,240,0)
 
    
    return C1 + H2 + H3 + H4
  

def Eth():
    
    '''
    takes in the eth group and creates the coordinates for group 
    
    the first carbon is at the coordinate 0,0,0 as standard 
    confusingly is not actually ethyl but vinyl
    '''
    
    
    C1 = (0,0,0,'C')
    C2 = (0,-Bond_Length.C_Cd,0,'C')
    H3 = (math.cos(math.radians(30))*Bond_Length.C_H,
                   -Bond_Length.C_Cd - math.sin(math.radians(30)*Bond_Length.C_H),
                   0,'H')
    H4 = (-math.cos(math.radians(30))*Bond_Length.C_H,
                    -Bond_Length.C_Cd - math.sin(math.radians(30)*Bond_Length.C_H),
                    0,'H')
    
    H5 = (-math.cos(math.radians(30))*Bond_Length.C_H,
          math.sin(math.radians(30))*Bond_Length.C_H,
          0,'H')
    
    return [C1,C2,H3,H4,H5]
    
    
    


def Backbone():
        
        '''
        creates the base Si 0 back bone 
        these will have the numbers 1 and 2 associated with them 
        
        they will be outputed as a list 
        '''
        
        Si1 = (0,0,0,'Si') 
        O2 = (round(math.cos(math.radians(35.25))*Bond_Length.Si_O,4),
              round(math.sin(math.radians(35.25))*Bond_Length.Si_O,4),
              0, 'O')
        
        return [Si1, O2] 


####---------------------------------------------------------------------------#####
####---------------------------------classes-----------------------------------#####
####---------------------------------------------------------------------------#####


    
class Polymer: 
    
    @property 
    def midpoint(self):
        coordinates = self.coordinates
        # print(len(self.coordinates))
        
        xlst = [values[0] for values in coordinates]
        ylst = [values[1] for values in coordinates]
        zlst = [values[2] for values in coordinates]
        
        xlo,xhi = min(xlst),max(xlst)
        ylo,yhi = min(ylst),max(ylst)
        zlo,zhi = min(zlst),max(zlst)    
                
        xmid = (xhi-xlo)/2 + xlo
        ymid = (yhi-ylo)/2 + ylo
        zmid = (zhi-zlo)/2 + zlo
        
        return np.array([xmid,ymid,zmid])
        
        
    @property 
    def radius(self):
        coordinates = self.coordinates
        
        xlst = [values[0] for values in coordinates]
        ylst = [values[1] for values in coordinates]
        zlst = [values[2] for values in coordinates]
        
        xlo,xhi = min(xlst),max(xlst)
        ylo,yhi = min(ylst),max(ylst)
        zlo,zhi = min(zlst),max(zlst)  
                
        xmid = (xhi-xlo)/2 
        ymid = (yhi-ylo)/2 
        zmid = (zhi-zlo)/2
        
        return (xmid**2 +ymid**2 +zmid**2)**.5
    
    @property 
    def box(self): 
        """
        Returns
        -------
        np.array showing np.array([[xmin,xmax],[ymin,yma...],...]

        """
        
        coordinates = self.coordinates
        
        xlst = [values[0] for values in coordinates]
        ylst = [values[1] for values in coordinates]
        zlst = [values[2] for values in coordinates]
        
        xlo,xhi = min(xlst),max(xlst)
        ylo,yhi = min(ylst),max(ylst)
        zlo,zhi = min(zlst),max(zlst) 
        
        return np.array([[xlo,xhi],[ylo,yhi],[zlo,zhi]]).T
    
    
    def __init__(self,coord_in):
        """

        Parameters
        ----------
        coord_in : list of tuples
            takes in a list of tuples and turns them into a polymer type object 
            example [(x,y,z,'atom_type'),(1,3,8,'C')]

        Returns
        -------
        returns the polymer class type object 

        """
        
        self.coordinates = coord_in
        
    def Rotate(self,psi,theta,phi):
        """

        Parameters
        ----------
        psi : int/float
            Rotation in degrees around the x axis 
        theta : int/float
            Rotation in degrees around the y axis
        phi : int/float 
            Rotation in degrees around the z axis

        Returns
        -------
        None. Changes the coordinate attribute values 

        """
        
        x = Rotation_func(self.coordinates,psi,theta,phi)
        self.coordinates = x 
        
    def Translate(self,trans):
        """

        Parameters
        ----------
        trans : tuple or list
            takes in a transformation vector of [x,y,z] or (x,y,z)

        Returns
        -------
        Updates the coordinate attribute

        """
        
        x = Translation_func(self.coordinates,trans)
        self.coordinates = x 
        
    def Mirror(self,plane):
        
    
        """

        Parameters
        ----------
        plane : str
            takes in str of 'x' 'y' or 'z' and reflects the coordinates over that plane

        Returns
        -------
        Updates the coordinate attribute

        """
        
        x = Mirror_func(self.coordinates,plane)
        self.coordinates = x   
    def Duplicate(self):
        """

        Returns
        -------
        TYPE : Polymer type obejct
            returns a polymer type object with a copied list

        """
        
        return self.__class__(self.coordinates.copy())
    
    def to_file(self,name,location,file_type='lammps'):
        """
        Parameters
        ----------
        name : str
            file name of output
        location : str
            path to the directory
        file_type : str, optional
            The default is 'lammps'.
            'xyz' is another possible input

        Returns
        -------
        creates a file and prints a string

        """
        
        tot_file = location +'/' + name
        file_maker(self.coordinates,tot_file,file_type)
        
        print(f'file created : {name}\nin dir: \n{location}')
        
    def pbc_translate(self,shift,box = None):
        
        # d = {'Si':1,'O':2,'C':3,'H':4}
        # d_rev = {1:'Si',2:'O',3:'C',4:'H'}
        # print(self.box.T)
        if box == None:
            box = self.box.T
        else:
            box = box.T 
        self.Translate(shift)
        coordinates = self.coordinates
        coord_matrix = []
        atom_ids = []
        for coords in coordinates: 
            x,y,z = coords[0:3]
            atom_id = coords[3]
            
            coord_matrix.append([x,y,z])
            atom_ids.append(atom_id)
            
        # print(coord_matrix,atom_ids)
        coord_matrix = np.array(coord_matrix).T
        
        # print(box)
        trans_matrix = []
        
        i = 0 
        while i < 3: 
            min_v = box[i][0]
            max_v = box[i][1]
            
            d = max_v - min_v 
            
            min_trans = (coord_matrix[i] < min_v).astype(int) * d 
            max_trans = (coord_matrix[i] > max_v).astype(int) * -d 
            
            trans = min_trans + max_trans 
            
            trans_matrix.append(trans)
            
            i += 1 
            
        trans_matrix = np.array(trans_matrix)
        
        # print(trans_matrix)
        new_coords = trans_matrix + coord_matrix
        # print(coord_matrix)
        new_coords = new_coords.T 
        # print(new_coords)
        
        pbc_coords = []
        i = 0 
        while i < len(atom_ids):
            temp = (new_coords[i][0],new_coords[i][1],new_coords[i][2],atom_ids[i])
            pbc_coords.append(temp)
            i += 1
            
        
        self.coordinates = pbc_coords
            
            
        
        
        
    @classmethod
    def file_read(cls,file_location,file_name,atom_coord=[1,2,3,4],file_in = False):
        """

        Parameters
        ----------
        cls : Polymer subclass or Polymer
            end data type you wish the coordinates to be 
        file_location : str
            Directory where file is stored
        file_name : str
            name of the file 
        atom_coord : list or tuple
            [Si,O,C,H] the numbers corresponding to them 
        file_in : bool (default: False) 
            if the file to read is an input file

        Returns
        -------
        returns a cls type object 

        """
        file_full = file_location + '/' + file_name
        
        if file_in == True: 
            d = {}
            d[atom_coord[0]] = 'Si'
            d[atom_coord[1]] = 'O'
            d[atom_coord[2]] = 'C'
            d[atom_coord[3]] = 'H'
            
            with open(file_full,'r') as f: 
                stuff = f.read() 
            
            stuff = stuff.split('\n')
            s_len = len(stuff)
            info_start = 0
            for i in range(len(stuff)): 
                if stuff[i][0:5] == "Atoms":
                    info_start = i 
                    break 
            stuff = stuff[info_start+2:s_len-1]
            # print(stuff)
            
            coords = []
            
            for values in stuff: 
                # print(values)
                con = values.split()
                # print(con)
                tup = (float(con[3]),float(con[4]),float(con[5]),d[int(con[1])])
                coords.append(tup)
                
            return cls(coords)
                
            
            
        
        d = {}
        d[atom_coord[0]] = 'Si'
        d[atom_coord[1]] = 'O'
        d[atom_coord[2]] = 'C'
        d[atom_coord[3]] = 'H'
        
        

        with open(file_full) as f: 
            stuff = f.read() 
            
           
        values = stuff.split('ITEM: TIMESTEP')
        # print(len(values))
        values = values[-1]
        values = values.split('\n')
        # print(values[0:10])
        # print(values[9])
        box = values[5:8]
        # print(box)
        
        temp_list = []
        for things in box :
            temp = things.split(' ')
            temp = [float(val) for val in temp]
            temp_list.append(temp)
            
        box = np.array([temp_list[0],temp_list[1],temp_list[2]])
        # print(box.T)
        ran = box.T[1] - box.T[0]
        # print(ran)
        
        coordinates = []
        raw_coord = values[9:-1]
        # print(raw_coord[-20:])
        
        for coord in raw_coord: 
            temp = coord.split(' ')
            # print(temp)
            xyz = np.array([float(temp[2]),float(temp[3]),float(temp[4])])
            xyz = xyz * ran 
            xyz = xyz + box.T[0]
            
            tup = (xyz[0],xyz[1],xyz[2],d[int(temp[1])])
            coordinates.append(tup)
            
        return cls(coordinates)
    
    def __add__(self,other):
        if issubclass(type(other),Polymer): 
            return self.__class__(self.coordinates + other.coordinates)
        else: 
            return self.__class__(self.coordinates + other)
        
    def __iadd__(self,other):
        if issubclass(type(other),Polymer): 
            self.coordinates += other.coordinates
            return self
        else: 
            self.coordinates += other
            return self
        
    def __str__(self):
        d = {'Si':0,'O':0,'C':0,'H':0}
        for values in self.coordinates: 
            d[values[3]] += 1
        total_atoms = len(self.coordinates)
        return f'total atoms : {total_atoms}' + '\n' + f'Si : {d["Si"]}' + '\n' + f'O : {d["O"]}' + '\n' + f'C : {d["C"]}' + '\n' + f'H : {d["H"]}'+'\n'
        



class linear_Polymer(Polymer):
    
    def __init__(self,coord_in,expansion_axis= np.array([1,0,0])):
        self.expansion_axis = expansion_axis
        self.coordinates = coord_in
                
        
    def Si_ORandCreator(length,side_group):
        
        
    
        '''
        Parameters
        ----------
        length : int 
            number of Si and O in the chains
        side_group : dict
        shows the side groups and their prevelence 
        {'Phenyl':2,'Methyl':2,'Vinyl':2,'Hydride':1}
        currently only 4 side groups and it will take the numbers and add them 
        then assign them the polymer side groups using random num generation            

        Returns
        -------
        None.

        '''
        
        ####side groups and main segments####
        Benzene_base = Polymer(Benzene())
        Benzene_base.Rotate(0,0,120)
        Ethyl_base = Polymer(Eth())
        Methyl_base = Polymer(Methal())
        
        BackBone_base = Polymer(Backbone())
        
        ####rearranged left and right base side groups####
        BenzeneR = Benzene_base.Duplicate()
        BenzeneR.Translate([0,-1.39,0])
        BenzeneR.Rotate(-90+35.25,60,0)
        BenzeneR.Translate([0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_C,math.cos(math.radians(35.25))*Bond_Length.Si_C])
        
        BenzeneL = Benzene_base.Duplicate()
        BenzeneL.Translate([0,-1.39,0])
        BenzeneL.Rotate(90-35.25,60,0)
        BenzeneL.Translate([0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_C,-math.cos(math.radians(35.25))*Bond_Length.Si_C])
        
        EthylL = Ethyl_base.Duplicate()
        EthylL.Rotate(0,90,0)
        EthylL.Translate([0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_C,-math.cos(math.radians(35.25))*Bond_Length.Si_C])
        
        EthylR = Ethyl_base.Duplicate()
        EthylR.Rotate( 0,-90,0)
        EthylR.Translate([0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_C,1 *math.cos(math.radians(35.25))*Bond_Length.Si_C])
        
        MethylL = Methyl_base.Duplicate()
        MethylL.Rotate(90-35.25, 0, 0)
        MethylL.Translate( [0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_C,-math.cos(math.radians(35.25))*Bond_Length.Si_C])
        
        MethylR = Methyl_base.Duplicate()
        MethylR.Rotate(-90+35.25, 0, 0)
        MethylR.Translate([0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_C,math.cos(math.radians(35.25))*Bond_Length.Si_C])
        
        HR = Polymer([(0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_H,math.cos(math.radians(35.25))*Bond_Length.Si_H,'H')])
        
        HL = Polymer([(0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_H,-math.cos(math.radians(35.25))*Bond_Length.Si_H,'H')])
        
        
        Components = [HR,HL,MethylL,MethylR,EthylL,EthylR,BenzeneL,BenzeneR,BackBone_base]
        
        
# =============================================================================
#         random num with associated with the side group 
# =============================================================================

        counter = 0 
        type_list = []
        for k, values in side_group.items() :
            # print(counter)
            counter += values 
            
            while values > 0:
                type_list.append(k)
                
                values -= 1 
                
        # print(counter)
                
        # print(type_list)
        length_list = range(counter)
        
        d = dict(zip(length_list,type_list))
        # print(d)
        
# =============================================================================
#     creating polymer first unit
# =============================================================================
            
        Base = BackBone_base.Duplicate()
        
        build_trans = [ 2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0]
        
        Lg = d[randint(0,counter-1)]
        Rg = d[randint(0,counter-1)]
        
        # print(Lg,Rg)
        
        #group sides 
        #{'Phenyl':2,'Methyl':2,'Vinyl':2,'Hydride':1}
        
        if Lg == 'Phenyl':
            Base += BenzeneL
        elif Lg == 'Methyl':
            Base += MethylL
        elif Lg == 'Vinyl':
            Base += EthylL
        else:
            Base += HL
            
        if Rg == 'Phenyl':
            Base += BenzeneR
        elif Rg == 'Methyl':
            Base += MethylR
        elif Rg == 'Vinyl':
            Base += EthylR
        else:
            Base += HR
        
        for i in range(length-1):
            
            for coms in Components:
                coms.Translate(build_trans)
                
            Base += BackBone_base
            
            Lg = d[randint(0,counter-1)]
            Rg = d[randint(0,counter-1)]
            
            # print(Lg,Rg)
            
            if Lg == 'Phenyl':
                Base += BenzeneL
            elif Lg == 'Methyl':
                Base += MethylL
            elif Lg == 'Vinyl':
                Base += EthylL
            elif Lg == 'Hydride':
                Base += HL
                
            if Rg == 'Phenyl':
                Base += BenzeneR
            elif Rg == 'Methyl':
                Base += MethylR
            elif Rg == 'Vinyl':
                Base += EthylR
            elif Rg == 'Hydride':
                Base += HR
                
                
        #adding two end groups 
        end_group1 = Polymer(Methal()) 
        end_group1.Rotate( 0, 180, 35.25+90)
        end_group1.Translate([-Bond_Length.Si_C*math.cos(math.radians(35.25)),
                                              Bond_Length.Si_C*math.sin(math.radians(35.25)),
                                              0])
        
        end_group2= Polymer(Methal()) 
        end_group2.Rotate( 0, 180, 35.25-90)
        end_group2.Translate([Bond_Length.Si_O*math.cos(math.radians(35.25)) * 2 * (length - .5) + 
                                          Bond_Length.C_O*math.cos(math.radians(35.25)),
                                          Bond_Length.Si_O * math.sin(math.radians(35.25))-
                                          Bond_Length.C_O * math.sin(math.radians(35.25)),
                                          0])
            
        Base += end_group1
        Base += end_group2
            
        return linear_Polymer(Base.coordinates)
            
    def PDES(length): 
        
        length = length // 2 + length % 2 
        
        def eth_true(): 
            
            '''
            Return a polymer class object that is an ethyl group 
            -------
        
            '''
            
            methal1 = Polymer(Methal())
            methal2 = methal1.Duplicate()
            methal1.Rotate(0,180,0)
            methal1.Rotate(-90,0,0)
            methal1.Translate([0,0,Bond_Length.C_C])
            methal2.Rotate(90,0,0) 
            
            z = methal1+methal2
            del z.coordinates[5]
            z.Rotate(180,0,0)
            z.Rotate(0,0,90)
            z.Rotate(-54.75,0,0)
            z.Rotate(0,180,0)
            return z
        
        ethyl_groupR= eth_true()
        ethyl_groupRR = ethyl_groupR.Duplicate()
        ethyl_groupRR.Rotate(60,180,0)
        ethyl_groupRR.Translate([2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                              0,
                                              0])
        R = ethyl_groupR + ethyl_groupRR
        R.Rotate(-90,0,0)
        R.Translate([0,-math.sin(math.radians(35.25))*Bond_Length.Si_C,math.cos(math.radians(35.25))*Bond_Length.Si_C])
        L = R.Duplicate()
        L.Mirror('z')
        T = R+L
        
        BackBone_base = Polymer(Backbone())
        base1 = BackBone_base.Duplicate()
        base2 = BackBone_base.Duplicate() 
        base2.Translate([2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,0,0])
        
        T += base1 + base2 
        m_T = T.Duplicate() 
        
        build_trans = [ 4 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0]
        
        for i in range(length-1): 
            m_T.Translate(build_trans)
            T+= m_T
            
        end_group1 = Polymer(Methal()) 
        end_group1.Rotate( 0, 180, 35.25+90)
        end_group1.Translate([-Bond_Length.Si_C*math.cos(math.radians(35.25)),
                                              Bond_Length.Si_C*math.sin(math.radians(35.25)),
                                              0])
        
        end_group2= Polymer(Methal()) 
        end_group2.Rotate( 0, 180, 35.25-90)
        end_group2.Translate([Bond_Length.Si_O*math.cos(math.radians(35.25)) * 2 * (length*2 - .5) + 
                                          Bond_Length.C_O*math.cos(math.radians(35.25)),
                                          Bond_Length.Si_O * math.sin(math.radians(35.25))-
                                          Bond_Length.C_O * math.sin(math.radians(35.25)),
                                          0])
            
        T += end_group1
        T += end_group2
        
        # T.to_file('PDES_expanded.in',r'G:\My Drive\Research\reaxff\lammps_in\PDES')
        
        return linear_Polymer(T.coordinates)
        


            
    def Rand_Rotation(self):
        """
        Returns
        -------
        None. - changes the self.coordinates of a linear polymer to 
        so that the polymer is rotated around some random coordinates

        """
        
        mid_trans = np.array(self.midpoint)
        
        psi, phi, theta = randint(0,360), randint(0,360), randint(0,360)
        
        self.Translate(-mid_trans)
        self.Rotate(psi,theta,phi)
        self.Translate(mid_trans)
        
    
    def __add__(self,other) :
        
        self.Rand_Rotation()
        other.Rand_Rotation() 
        
        self_exp = self.box * self.expansion_axis
        other_exp = other.box * self.expansion_axis
        
        # print(self.box)
        # print(other.box)
        
        self_exp = self_exp[1]
        other_exp = other_exp[0]
        
        trans = self_exp - other_exp 
        other.Translate(trans)
        
        expansion = self.expansion_axis
        print(expansion)
        exp_matrix = np.array([[0,1,0],[0,0,1],[1,0,0]]).T
        expansion = exp_matrix @ expansion
        
        coords = self.coordinates + other.coordinates 
        
        return linear_Polymer(coords,expansion)
        
        
        
        
class amorphous_box(Polymer):
    
    def __init__(self,coord_in,expansion_axis=np.array([1,0,0])):
        self.expansion_axis = expansion_axis
        self.coordinates = coord_in
        
    def rand_expansion(self,shift,direct=None,buffer=.07,dup_num=1):
        """
        
        Parameters
        ----------
        shift : tuple 
            direction of pbc shift and how much to shift by 
        direct : tuple or default = None 
            if None then will expand in the x direction 
            if tuple will expand in that direction 
        buffer : float default = .07
            distance between the two boxes as a fraction of the box  

        Returns
        -------
        adds returns a amorphous type object with the coordinates doubled and pbc shifted

        """        
        if direct == None:
            expansion_axis = self.expansion_axis
        else: 
            expansion_axis = np.array(direct) 
            
        box = self.box
        d = box[1]- box[0]
        # print(d)
        d = d*expansion_axis*(1+buffer)
        # print(d)
        # print(d)
        
        T = self.Duplicate() 
        
        for i in range(dup_num): 
            exp = self.Duplicate()
            exp.Translate(d*(i+1))
            exp.pbc_translate(shift)
            T += exp 
            
        
        
        return T
        
        
            
        

        
        
        


####---------------------------------------------------------------------------#####
####---------------------------------tests-------------------------------------#####
####---------------------------------------------------------------------------#####




          
# print(Rotation([(1,2,3,'C'),(1,8,9,'Si')],90,90,90))
# print(Translation([(1,2,3,'C'),(1,8,9,'Si')],(1,2,3)))
# print(Mirror([(1,2,3,'C'),(1,8,9,'Si')],'x'))

# print(Benzene())
# print(Methal())
# print(Eth())
      

# Benzene1 = Polymer(Benzene())
# Benzene1.Rotate(0,0,120)

# Benzene2.Translate((3,3,3))
# BenzeneT = Benzene1 + Benzene2

# Benzenef = Benzene1.Duplicate()
# Benzenef += Benzene2
# Benzene1.Translate([0,-1.39,0])
# # print(Benzene1.coordinates)
# Benzene2 = Benzene1.Duplicate()
# Benzene3 = Benzene1.Duplicate()
# Benzene1.to_file('Benzene_test1.in','G:/My Drive/Research/reaxff/lammps_in')
# Benzene1.Rotate(-90+35.25,60,0)
# Benzene1.to_file('Benzene_test2.in','G:/My Drive/Research/reaxff/lammps_in')
# z = Benzene1 
# Benzene2.Rotate(0,60,0)
# Benzene3.Rotate(-90+35.25,0,0)
# z += Benzene3
# z.to_file('Benzene_test3.in','G:/My Drive/Research/reaxff/lammps_in')



# file_maker(coords,'G:/My Drive/Research/reaxff/lammps_in/Benzene_test.in')
        


# [print(coords) for coords in Benzenef.coordinates]
# print([Benzenef,Benzene2])


# x = linear_Polymer.Si_ORandCreator(100, {'Phenyl':0,'Methyl':0,'Vinyl':2,'Hydride':1})
# x.Rand_Rotation()
# print(x.midpoint)
# x.Translate((20,20,20))
# print(x.midpoint)
# trans = np.array([1,1,1])
# x.Translate(trans)
# x.to_file('lin_test_2-9-1.in','G:/My Drive/Research/reaxff/lammps_in' )
# print(x)

# x = linear_Polymer.Si_ORandCreator(40, {'Phenyl':1,'Methyl':1,'Vinyl':0,'Hydride':0})
# y = linear_Polymer.Si_ORandCreator(40, {'Phenyl':0,'Methyl':1,'Vinyl':0,'Hydride':1})
# n = linear_Polymer.Si_ORandCreator(40, {'Phenyl':0,'Methyl':1,'Vinyl':0,'Hydride':1})
# m = linear_Polymer.Si_ORandCreator(40, {'Phenyl':1,'Methyl':1,'Vinyl':0,'Hydride':0})
# p = linear_Polymer.Si_ORandCreator(40, {'Phenyl':1,'Methyl':1,'Vinyl':0,'Hydride':0})
# print(x.box)
# z = x + y + n + m + p
# print(z.box)
# z.to_file('lin_test_2-24-5.in','G:/My Drive/Research/reaxff/lammps_in')
# y = x.midpoint
# x.Rand_Rotation()
# print(x.box)
# print(x.box[0][1])
# print(x.midpoint-y)

# x = amorphous_box.file_read('G:/My Drive/Research/reaxff/mc/PSO_t8_5','loop_heat.1',[3,4,2,1])
# print(x)
# x.to_file('end_group_test1.in','G:/My Drive/Research/reaxff/lammps_in')
# z = x.rand_expansion((-20,20,30))
# print(z)
# z.to_file('expand_test1.in','G:/My Drive/Research/reaxff/lammps_in')
# print(x)
# x.to_file('read_test.in','G:/My Drive/Research/reaxff/lammps_in')

# x = linear_Polymer.Si_ORandCreator(100, {'Phenyl':0,'Methyl':0,'Vinyl':2,'Hydride':1})
# x.Rand_Rotation()
# x.to_file('pbc_test1.in','G:/My Drive/Research/reaxff/lammps_in')
# x.pbc_translate((0,0,20))
# x.to_file('pbc_test2.in','G:/My Drive/Research/reaxff/lammps_in')
# x.pbc_translate((0,0,20))
# x.to_file('pbc_test3.in','G:/My Drive/Research/reaxff/lammps_in')
# x.pbc_translate((0,0,20))
# x.to_file('pbc_test4.in','G:/My Drive/Research/reaxff/lammps_in')
# x.pbc_translate((0,0,20))
# x.to_file('pbc_test5.in','G:/My Drive/Research/reaxff/lammps_in')
# print(x.coordinates)
# print(x)

# x = linear_Polymer.Si_ORandCreator(100, {'Phenyl':0,'Methyl':0,'Vinyl':2,'Hydride':1})
# x.Rand_Rotation()
# bx = x.box
# for value in range(100):
#     x.to_file(f'pbc_movie{value}.in',r'G:\My Drive\Research\reaxff\lammps_in\pbc_test_movie')
#     x.pbc_translate((5,-6,5),bx)


# x = amorphous_box.file_read('G:/My Drive/Research/reaxff/mc/3d_relax_t7_1','loop_heat.25',[1,2,4,3])
# x.pbc_translate((30,0,-10))
# x.to_file('pso_strength.in','G:/My Drive/Research/reaxff/lammps_in')


# x = amorphous_box.file_read(r'G:\My Drive\Research\reaxff\PDMS','PDMS_t2_2.rdx',[1,2,3,4])
# x = x.rand_expansion((0,0,0),(1,0,0),buffer=.04)
# x = x.rand_expansion((0,0,0),(0,1,0))
# x = x.rand_expansion((0,0,0),(0,0,1))
# # x = x.rand_expansion((0,0,0),(0,1,0),buffer=.04)
# x.to_file('PDMS_expanded.in','G:/My Drive/Research/reaxff/lammps_in/PDMS')

# x = amorphous_box.file_read(r'G:\My Drive\Research\reaxff\lammps_in\PPSQ',r'PPSQ.in',[1,2,3,4],file_in=True)
# print(x.box)
# x = x.rand_expansion((0,0,0),(1,0,0),buffer=.13)
# x = x.rand_expansion((0,0,0),(0,1,0),buffer=.17)
# x = x.rand_expansion((0,0,0),(0,0,1),buffer=.13)
# x = x.rand_expansion((0,0,0),(1,0,0),buffer=.13)
# x = x.rand_expansion((0,0,0),(0,1,0),buffer=.13)
# x = x.rand_expansion((0,0,0),(0,0,1),buffer=.13)
# x = x.rand_expansion((0,0,0),(1,0,0),buffer=.13)
# x = x.rand_expansion((0,0,0),(0,1,0),buffer=.13)
# x = x.rand_expansion((0,0,0),(0,0,1),buffer=.13)
# print(x)
# x.to_file('PPSQ_expanded.in',r'G:\My Drive\Research\reaxff\lammps_in\PPSQ')


# for i in range(20): 
#     x1 = linear_Polymer.PDES(100)
#     x2 = linear_Polymer.PDES(100)
#     x3 = linear_Polymer.PDES(100)
#     x4 = linear_Polymer.PDES(100)
#     x5 = linear_Polymer.PDES(100)
#     x6 = linear_Polymer.PDES(100)
#     x7 = linear_Polymer.PDES(100)
#     x8 = linear_Polymer.PDES(100)
    
#     z = x1+x2+x3+x4+x5+x6+x7+x8
#     z.to_file(f'PDES_randex.{i}',r'G:\My Drive\Research\reaxff\lammps_in\PDES\test')
