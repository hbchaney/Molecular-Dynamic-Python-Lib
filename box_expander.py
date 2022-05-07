# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 14:33:16 2021

@author: hbcha
"""

'''
unit cell expander 

needs: 
    enough atoms to observe macro strucutres 
    ovserve 3d structures 
     
constraint: 
    unchanged proportions 
    
routes:
    
creating a duplicate cell next to the existing cell 
     -would keep the same proportions 
     -limits control over # of atoms 
     -might just further the existing 2d structure 
     
creating a duplicate cell that is then rotated  <--------------------
    -same proportions 
    -same limit on # of atoms (factors of whole integers)
    -might cause significant instability 
    -would have some kind of change on the structure 
    
expanding the unit cell by some factor less than 1
    -diffcult to keep the at fractions the same 
    -more control over the number of atoms 
    -might further the 2d structure
    
expanding by factor less than 1 and rotating 
    -difficult to implement a repressentative region would need to be chosen 
    -more control over # of atoms 
    -would change up the structure
    -might introduce too much instability 
    
alternative routes:
    
micro sampling
    - find a recatngular region where the atom fraction is repressentative of the whole roughly 
    than create duplicates that are then stacked randomly throughout the box 
    -not truly from the pyrolosis route and might not be repressentative 
    -would definately show the 3d behavior 
        
'''



'''
inputs: 
    structure file in the lammps output format 
    
outputs:
    an xyz file to be inputed into ovito 

'''

###############--------------------------------------------------------------------################

'''
imports
'''

# import matplotlib.pyplot as plt 
import math 
import numpy as np 
from decimal import Decimal 



################-----------------------------------------------------------------######################
######### New Functions #############################################################################
################-----------------------------------------------------------------####################

####### rdx coord finder ###########-----------------------------------------------------------------

def rdx_reader(filename,path):
    '''
    takes in a file from a specified path 
    outputs the coordinate data in the standard polymer format 
    as well as some box data 
    [x,y,z,ele,id] ex [1.01,1.05,1.06,'c',56]

    '''
    #input 0 if the path is in the cwd 
    #input the wd of the filename if not in the cwd 
    
    if path != 0 :
        with open(path + '/' + filename) as f:
            contents = f.readlines()
    else: 
        with open(filename) as f: 
            contents = f.readlines()
            
    # print(contents[0:10])
    # print(contents[0][6])
    
    ts = []
    counter = 0
    
    for lines in contents: 
        if len(lines) > 7:
            if lines[6] == 'T':
                ts.append(counter)
        
        counter += 1 
        
    box_info = contents[5:8]
    
    x = box_info[0].split()
    y = box_info[1].split()
    z = box_info[2].split()
    
    x[0], x[1] = float(Decimal(x[0])), float(Decimal(x[1]))
    y[0], y[1] = float(Decimal(y[0])), float(Decimal(y[1]))
    z[0], z[1] = float(Decimal(z[0])), float(Decimal(z[1]))
    
    # print(f'{x},{y},{z}')
    
    # print(box_info) 

    dx,dy,dz = x[1]-x[0],y[1]-y[0],z[1]-z[0]
    # print(f'{dx} {dy} {dz}')
    # print(ts)
    # print(contents[3858:3863])
    contents = contents[9:ts[1]]
    # print(contents[-1])
    
    type_dict = {'1':'H','2':'C','3':'Si','4':'O'}
    
    atom_coords = []
    
    for values in contents: 
        sep = values.split(' ')
        x,y,z = float(sep[2])*dx,float(sep[3])*dy,float(sep[4])*dz
        ele = type_dict[sep[1]]
        atm_id = int(sep[0])
        
        atom = [x,y,z,ele,atm_id]   
        atom_coords.append(atom)
        
    return [[dx,dy,dz],atom_coords]
        
####### copy flipper ########-------------------------------------------------------------------

def flip_copy(axis,trans,direct,coords,box):
    '''
    takes a list of coords and copys and extends them across an axis 
    it then rotates the coordinates by 180 degrees along the flip axis 
    returns both set of coords
    '''
    
    coord_copy = coords.copy()
    
    dx,dy,dz = box
    dx2,dy2,dz2 = dx/2,dy/2,dz/2
    buffer = 2
    
    
    #moving in the direction of axis 
    if axis == 'x':
        ax_move = [dx+buffer,0,0]
        ax_rot = [180,0,0]
    elif axis == 'y':
        ax_move = [0,dy+buffer,0]
        ax_rot = [0,180,0]
    else: 
        ax_move = [0,0,dz+buffer]
        ax_rot = [0,0,180]
        
    coord_copy = Translation(coord_copy,[-dx2,-dy2,-dz2])
    coord_copy = Rotation(coord_copy,ax_rot[0],ax_rot[1],ax_rot[2])
    coord_copy = Translation(coord_copy,[dx2,dy2,dz2])
    coord_copy = Translation(coord_copy,ax_move)
    
    coord_good = []
    coord_bad = []
        
    if direct == 'x':
        shift = [trans,0,0]
        coord_copy = Translation(coord_copy,shift)
        for coord in coord_copy: 
            if coord[0] > dx :
                coord_bad.append(coord)
            else: 
                coord_good.append(coord)
                
        coord_bad = Translation(coord_bad,[-dx,0,0])
    elif direct == 'y':
        shift = [0,trans,0]
        coord_copy = Translation(coord_copy,shift)
        for coord in coord_copy: 
            if coord[1] > dy :
                coord_bad.append(coord)
            else: 
                coord_good.append(coord)
                
        coord_bad = Translation(coord_bad,[0,-dy,0])
    else:
        shift = [0,0,trans]
        coord_copy = Translation(coord_copy,shift)
        for coord in coord_copy: 
            if coord[2] > dz :
                coord_bad.append(coord)
            else: 
                coord_good.append(coord)
                
        coord_bad = Translation(coord_bad,[0,0,-dz])
    

    
    
    print(len(coord_good))
    print(len(coord_bad))
    print(len(coords))
    coord_good = coord_add(coord_good,coords)
    coord_bad = coord_add(coord_bad,coord_good)
    
    # print((coords + coord_good + coord_bad)[-20:])
    return coords + coord_good + coord_bad
    

###### previously made functions #######--------------------------------------------------------


#### xyz name creator #### -------------------------------------------------------------

def xyz_name(filename): 
    
    '''creates the xyz file with the correct file name
    also sets up some important preliminary stuff 
    like the header lines'''
    
    with open(filename,'w') as f:
        f.write('number\n')
        f.write('notes\n')
    
    print(f'file created {filename}') 
    
#### xyz appender ####-----------------------------------------------------------------

def xyz_appender(coord,filename):

    
    '''
    takes in the list of coordinates and formats them into a lines to add to the xyz file 
    and formats then adds them 
    
    ''' 
    xyz_name(filename)
    
    
    with open(filename,'a') as f: 
        for coords in coord: 
            f.write(f'{coords[3]}  {coords[0]}  {coords[1]}  {coords[2]} \n')
            
    atom_id_list = []
    for values in coord: 
        atom_id = values[4]
        atom_id_list.append(atom_id)
        
    max_id = len(coord)
    last_num = str(max_id) + '\n'
            
    with open(filename,'r') as f: 
        lines = f.readlines() 
        
    lines[0] = last_num 
    
    with open(filename, 'w') as f: 
        f.writelines(lines) 

#### Translation function ####----------------------------------------------------------

def Translation(coord_raw,trans): 
    
    '''
    takes in a list of coordinates and and specific translation and outputs 
    the a list translated coordiantes 
    
    can also output individual coordinates if needed 
    
    trans input: 
        [x,y,z]
    '''
    
    trans = np.array(trans)
    
    #check to see if the it is one coordinate or a list of coordinates 
    if type(coord_raw[0]) != list: 
        
        #first parsing to get only the coordinate info 
        coord_array = np.array(coord_raw[0:3]) 
        coord_new_array = trans + coord_array
        #making the location back into a list 
        coord_new_unrounded = coord_new_array.tolist()
        
        #rounding again 
        coord_new = []
        for values in coord_new_unrounded: 
            v_round = round(values, 4)
            coord_new.append(v_round)
        
        #adding back the extra info taken off 
        coord_new.append(coord_raw[3])
        coord_new.append(coord_raw[4])
        
        
        
        return coord_new
    
    else: 
        
        coord_list = []
        
        for coords in coord_raw: 
            
            coords_array = np.array(coords[0:3])
            coords_new_array = trans + coords_array
            coords_new_unrounded = coords_new_array.tolist()
            
            #rounding again 
            coords_new = []
            for values in coords_new_unrounded: 
                v_round = round(values, 4)
                coords_new.append(v_round)
            
            coords_new.append(coords[3])
            coords_new.append(coords[4])
            
            coord_list.append(coords_new)
            
        return coord_list

#### Mirror function ####---------------------------------------------------------------

def Mirror(coord_raw,plane): 
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
    
    if type(coord_raw[0]) != list: 
        
        coord_array = np.array(coord_raw[0:3]) 
        coord_new_array = plane * coord_array
        coord_new = coord_new_array.tolist()
        coord_new.append(coord_raw[3])
        coord_new.append(coord_raw[4])
        
        return coord_new
    
    else: 
        
        coord_list = []
        
        for coords in coord_raw: 
            
            coords_array = np.array(coords[0:3])
            coords_new_array = plane * coords_array
            coords_new = coords_new_array.tolist()
            coords_new.append(coords[3])
            coords_new.append(coords[4])
            
            coord_list.append(coords_new)
            
        return coord_list

#### Rotation function ####------------------------------------------------------------

def Rotation(coord_raw,psi,theta,phi): 
    
    '''
    takes in a list of coordinates and rotates them about a specific axis(plural)
    can also take in a singular coordinate and output a singular value
    
    all rotations are clockwise?
    all rotation inputs in degrees 
    
    psi = rotation about the x axis 
    theta = rotation about the y axis 
    phi = rotation about the z axis 
    
    '''
    #creating the matricies of rotation 
    psi_matrix = np.matrix([[1, 0, 0], 
                           [0, math.cos(math.radians(psi)), math.sin(math.radians(psi))],
                           [0, -1*math.sin(math.radians(psi)), math.cos(math.radians(psi))]])
    
    theta_matrix = np.matrix([[math.cos(math.radians(theta)), 0, math.sin(math.radians(theta))],
                              [0, 1, 0],
                              [-1*math.sin(math.radians(theta)), 0, math.cos(math.radians(theta))]])
    
    phi_matrix = np.matrix([[math.cos(math.radians(phi)), math.sin(math.radians(phi)), 0],
                            [-1 * math.sin(math.radians(phi)), math.cos(math.radians(phi)), 0], 
                            [0, 0, 1]])
                            
    #combining the matricies to make a general rotation matrix 
    combined_matrix = np.matmul(np.matmul(psi_matrix,theta_matrix),phi_matrix)
    # print(combined_matrix)
    
    if type(coord_raw[0]) != list: 
        
        coord_array = np.array(coord_raw[0:3])
        new_coord_array = np.array(np.matmul(combined_matrix,coord_array)).tolist()[0]
        new_coord_unrounded = new_coord_array
        
        new_coord = []
        for values in new_coord_unrounded: 
            v_round = round(values, 4)
            new_coord.append(v_round)
        
        
        new_coord.append(coord_raw[3])
        new_coord.append(coord_raw[4])
        
        return new_coord
        
    
    else: 
        
        coord_list = []
        
        for coords in coord_raw: 
            
            coords_array = np.array(coords[0:3])
            new_coords_array = np.array(np.matmul(combined_matrix,coords_array)).tolist()[0]
            new_coords_unrounded = new_coords_array
            
            new_coords = []
            for values in new_coords_unrounded: 
                v_round = round(values, 4)
                new_coords.append(v_round)
            
            
            new_coords.append(coords[3])
            new_coords.append(coords[4])
            coord_list.append(new_coords)
        
        return coord_list
            
#### adder #### ----------------------------------------------------------------------

def coord_add(new_coords,last_num) :
    '''
    takes in the last number and adds it to all the numbers on the coords 
    
    
    the number is always in the 4th position 
    
    
    '''
    
    new_coords_list = new_coords.copy()
    
    if type(last_num[0]) == list: 
        atom_id_list =[]
        for coords in last_num: 
            atom_id = coords[4]
            atom_id_list.append(atom_id)
        last_number = max(atom_id_list)
        
    else: 
        # last_number = last_num[-1] 
        id_list = []
        for values in last_num:
            id_list.append(last_num[3])
        
        last_number = max(id_list)
        print(last_number)

    
    new_coords_updated = [] 
    
    if type(new_coords_list[0]) == list:

        for values in new_coords_list: 
            total_info = values.copy() 
            total_info[4] = total_info[4] + last_number
            new_coords_updated.append(total_info) 
    else: 
        new_coords_list[4] = new_coords_list[4] + last_number
        new_coords_updated = new_coords_list
    return new_coords_updated 
    
    
#############################################################################################

def main(): 
    
    filename = 'loop_heat.1'
    path = 'C:/Users/hbcha/Desktop/Researching work/Lammps PBC/python_gas/stable_pyro_t7_6'
    atom_og = rdx_reader(filename,path)
    # print(atom_og[0:10])
    
    coords = flip_copy('x',18,'z',atom_og[1],atom_og[0])
    xyz_appender(coords,'box_expander1.xyz')
    
    
    
    pass


if __name__ == '__main__':
    main()



















