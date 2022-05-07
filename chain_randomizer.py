# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 10:11:11 2021

@author: Harrison
"""

'''
the purpose of this program is to take a list of assigned values and create a box of 
n size with each box having a random orentation 

'''
#### Imported functions ####-----------------------------------------------------------
import numpy as np
import math 
import decimal
import random

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
    
#### adder #### ----------------------------------------------------------------------

def coord_add(new_coords,old_coords) :
    '''
    takes in the last number and adds it to all the numbers on the coords 
    
    
    the number is always in the 4th position 
    
    
    '''
    
    new_coords_list = new_coords.copy()
    '''
    if type(last_num[0]) == list: 
        last_copy = last_num.copy()
        last_number = last_copy[-1][4]
        
    else: 
        last_number = last_num[-1] 
    '''
    id_list = [] 
    for values in old_coords :
        id_list.append(values[4])
    
    last_number = max(id_list)

    
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
          
    last_num = str(len(coord)) + '\n'
            
    with open(filename,'r') as f: 
        lines = f.readlines() 
        
    lines[0] = last_num 
    
    with open(filename, 'w') as f: 
        f.writelines(lines) 
        
        
#######################################################################################
###### New functions ##################################################################

##### coord multiplyer #####------------------------------------------------------------
def coord_mult(coords_box): 
    '''
    takes in the coords and the box in string form  and multiplies them by the correct scale 
    also formats the box to be more easily read later on will be used within the rdx_parse function
    
    the output will be the same as the rdx only correctly formatted
    '''
    
    coords_raw = coords_box[0]
    box_raw = coords_box[1]
    
    
    box =[]
    for xyz in box_raw:
        xyz = xyz.split(' ')
        xyz = [float(decimal.Decimal(value)) for value in xyz]
        box.append(xyz)
        
    xscale = box[0][1] - box[0][0]
    yscale = box[1][1] - box[1][0]
    zscale = box[2][1] - box[2][0]
    
    # print(xscale,yscale,zscale)
    
    # print(box)
    new_coords = []
    for coords in coords_raw: 
        templist = []
        x = coords[0] * xscale
        y = coords[1] * yscale 
        z = coords[2] * zscale
        templist = [x,y,z,coords[3],coords[4]]
        new_coords.append(templist)
        
    return [new_coords,box]
    

#### rdx parser #### -------------------------------------------------------------------
def rdx_parse(filename,timestep): 
    '''
    takes in the filename of an rdx file and parses the line information in a specific time step 
    the file format is atom id , type, x,y,z 
    coords are [x,y,z,type,atom id] <--- foramt we need the info in 
    
    there are also many time steps in the files and they need to be distinguished 
    needs to also be general for any number of atoms 
    
    returns a list of two values 
    [[x,y,z,type,atom id],[box dimensions in str form]]
    
    '''
    
    with open(filename) as f :
        lines = f.readlines() 
        
    time_breaks = [] #format will be [[ts, linenum]
        
    #print(lines[0:20])
    
    #splitting up the code into sections based on the timestep 
    #timestep will be the same as show (0 will be included) 
    
    iter_vari = 0
    for values in lines: 
        #print(values[0:8])
        
        if values[0:8] == 'ITEM: TI': 
            temp_ts = [lines[iter_vari+1],iter_vari+1]
            time_breaks.append(temp_ts)
            
        iter_vari += 1
        
    #print(time_breaks)
    timestep = str(timestep) + '\n'
    time_range = []
    
    iter_vari = 0 
    for times in time_breaks: 
        
        if times[0] == timestep and times != time_breaks[-1]: 
            time_range.append([time_breaks[iter_vari][1],time_breaks[iter_vari+1][1]])
        
        if times == time_breaks[-1] and times[0] == timestep: 
            time_range.append([time_breaks[iter_vari][1],len(lines)])
        
        iter_vari += 1
    
    time_range = time_range[0]
    #print(time_range)
    position_raw = lines[time_range[0]:time_range[1]]
    #print(position_raw[-1])
    # print(position_raw[0:10])
    end_position = len(position_raw)-1
    # print(position_raw[end_position])
    positions = position_raw[8:end_position]
    # print(positions[-1])
    # print(positions[0:10])
    
    position_range = position_raw[4:7]
    position_range = [values[:-1] for values in position_range]
    #print(position_range)
    
    refined_pos =[]
    #rearranging the positions 
    for values in positions:
        values = values[:-1]
        temp_list = values.split(' ')
        temp_list = [float(values) for values in temp_list]
        #print(temp_list)
        
        x, y, z = temp_list[2],temp_list[3],temp_list[4]
        # print(x,y,z)
        atom_id = int(temp_list[0])
        atom_type = 'unkwn'
        
        if temp_list[1] == 1:
            atom_type = 'Si'
        elif temp_list[1] == 2:
            atom_type = 'O'
        elif temp_list[1] == 3: 
            atom_type = 'C'
        else:
            atom_type = 'H'
            
        refined_pos.append([x,y,z,atom_type,atom_id])
        
    
    
    return [refined_pos,position_range]
            
        
#### general box rotation randomizer ####----------------------------------------------
def rand_rotator(coords_box):
    '''
    takes in the coordinates oh the chain and rotates them randomly 
    also takes in the raw box data in the form of scientific notation 
    
    takes in the values with this format [refined_pos,position_range]
    where the refined positions = [x,y,z,type,atom id]
    
    returns [[x,y,z,type,atom id],boxsize element]
    
    '''
    # print(coords_box)
    box = coords_box[1]
    coords = coords_box[0]
    
    # print(box)
    #print(coords[0:20])
    
    #figuring out the center point and adjusting the values to make them roatate around the center
    

    
    # print(coords_box[1])
    
    center = []
    for values in box: 
        center_coor = (abs(values[1] - values[0]))/2
        center.append(-(center_coor))
    
    #print(center)
    
    #translating the to the center 
    
    #print(coords[0:10])
    coords = Translation(coords,center)
    #print(coords[0:10])
    
    psi_rand = random.randrange(0,360,1)
    theta_rand = random.randrange(0,360,1)
    phi_rand = random.randrange(0,360,1)
    coords = Rotation(coords, psi_rand,theta_rand,phi_rand)
    # print(coords[0:10])
    return coords
    
#### 3d random orentation box creation ####--------------------------------------------
def rand_box(dim,coords_box):
    '''
    this takes advantage of the rand_rotator and other functions to 
    create a box of randomly oreiented chains 
    
    the logic is as follows it creates a small box than moves and creates another box 
    rinse and repeat 
    
    return coords
    '''
    coords = coords_box[0]
    
    
    add_fact = 5.0
    box = coords_box[1] 
    max_len = box[0][1] - box[0][0] + add_fact
    
    i = 0 
    j = 0 
    k = 0 
    
    coords_list_com = [] 
    coords_list = []
    
    while i < dim: 
        j = 0
        while j < dim: 
            k = 0
            while k < dim: 
                
                if len(coords_list_com) != 0: 
                    coords = rand_rotator(coords_box)
                    coords = coord_add(coords,coords_list_com[-1])
                else: 
                    coords = rand_rotator(coords_box)
                    
                x = max_len * i
                y = max_len * j 
                z = max_len * k 
                
                xyz = [x,y,z]
                
                coords = Translation(coords,xyz)
                
                coords_list_com.append(coords)
                
                k += 1
            
            j += 1 
        
        i += 1
    
    
    for blocks in coords_list_com:
        for coords in blocks: 
            coords_list.append(coords)
            
    return coords_list

        
#### main fucntion ####----------------------------------------------------------------

def main(): 
    
    input_location = "G:/Shared drives/1 Kathy Lu Group/Harrison/reaxff/highcarbon/"
    input_name = 'high_carbon_t2.rdx'
    
    x = rdx_parse(input_location+input_name,80000)
    x = coord_mult(x)
    x = rand_box(2,x)
    
    xyz_name = 'hc_2x2_box_t2.xyz'
    xyz_dir = 'G:/Shared drives/1 Kathy Lu Group/Harrison/reaxff/xyz_files/'
    
    xyz_appender(x,xyz_dir+xyz_name)
    
    
    
if __name__ == '__main__': 
    main()




















