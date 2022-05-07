# -*- coding: utf-8 -*-
"""
Created by: Harrison Chaney 6/30/2021 

this will create polysiloxone xyz file that could be viewed using ovito 
and later converted into a lammps file 

many of the functions used can also be used for other  work 

notes: 
    all coord inputs are in this format [x,y,z,El,atom#]
    for example a carbon atom at the position (1,2,3) 
    would be [x,y,z,C,2]
    bond_list = [bond #, bond type,[atoms linked together]]
    angle_list = [angle #, angle ,[two atoms bonded]]
    

"""
#### imported modules ####---------------------------------------------------------
import numpy as np 
import math 
import random 




#### Global Variables ####----------------------------------------------------------

#global variables if any 
atom_counter = 0 
bond_counter = 0 
dihedral_counter = 0 


#### creating a class to hold the atom info ####-------------------------------------------

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
            
    
    # print(psi_matrix) 
    
    
    
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
        last_number = last_num[-1] 

    
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
    
    
    
    
#### Benzene creator ####---------------------------------------------------------------

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
    
    C3 = [(.5 * Bond_Length.C_CPh)/(math.tan(math.radians(30))), .5 * Bond_Length.C_CPh, 0, 'C', 3]
    H8 = [(C3[0] + math.cos(math.radians(30))), C3[1] + math.sin(math.radians(30)), 0, 'H', 8] 
    
    
    C3, H8 = Rotation([C3,H8],0,0,0)
    
    C5, H10 = Rotation([C3,H8],0,0,60)
    C5[4], H10[4] = 5, 10
    
    C6, H11 = Rotation([C3,H8],0,0,120)
    C6[4], H11[4] = 6, 11 
    
    C4, H9  = Rotation([C3,H8],0,0,180)
    C4[4], H9[4] = 4, 9 
    
    C2, H7  = Rotation([C3,H8],0,0,240)
    C2[4], H7[4] = 2, 7 
    
    C1 = Rotation(C3,0,0,300)
    C1[4] = 1
    
    # H12 = Rotation(H8,0,0,300) #take out if needed 
    # H12[4] = 12
    
    return [C1,C2,C3,C4,C5,C6,H7,H8,H9,H10,H11] #take out H12

#### CH3 creator function ####---------------------------------------------------------

def Methal(): 
    
    '''
    creates a methal group by using an initial point and rotating around the main bonded axis 
    
    the carbon created exists at the point [0,0,0]
    
        
    '''
    C1 = [0,0,0,'C',1]
    
    H2 = [math.cos(math.radians(19.5))*Bond_Length.C_H,
          -1 * math.sin(math.radians(19.5)) * Bond_Length.C_H,
          0,
          'H',
          2]
    
    H3 = Rotation(H2,0,120,0)
    H3[4] = 3
    
    H4 = Rotation(H2,0,240,0)
    H4[4] = 4 
    
    return [C1,H2,H3,H4]


#### Ethane group creator ####---------------------------------------------------------

def Eth() :
    
    '''
    takes in the eth group and creates the coordinates for group 
    
    the first carbon is at the coordinate 0,0,0 as standard 
    '''
    
    
    C1 = [0,0,0,'C',1]
    C2 = [0,-Bond_Length.C_Cd,0,'C',2]
    H3 = [math.cos(math.radians(30))*Bond_Length.C_H,
                   -Bond_Length.C_Cd - math.sin(math.radians(30)*Bond_Length.C_H),
                   0,'H',3]
    H4 = [-math.cos(math.radians(30))*Bond_Length.C_H,
                    -Bond_Length.C_Cd - math.sin(math.radians(30)*Bond_Length.C_H),
                    0,'H',4]
    
    H5 = [-math.cos(math.radians(30))*Bond_Length.C_H,
          math.sin(math.radians(30))*Bond_Length.C_H,
          0,'H',5]
    
    return [C1,C2,H3,H4,H5]


#### back bone creator #### -----------------------------------------------------------

def Backbone(): 
    
    '''
    creates the base Si 0 back bone 
    these will have the numbers 1 and 2 associated with them 
    
    they will be outputed as a list 
    '''
    
    Si1 = [0,0,0,'Si', 1] 
    O2 = [round(math.cos(math.radians(35.25))*Bond_Length.Si_O,4),
          round(math.sin(math.radians(35.25))*Bond_Length.Si_O,4),
          0, 'O', 2]
    
    return [Si1, O2] 

#### Segementor functions #####---------------------------------------------------------

def SiPh(): 
    '''
    creates a SiPh segment with xyz information

    the output will always be assuming first segment
    '''
    
    #the SiPh_list_comp will store the corrected components 
    SiPh_list_comp = []
    SiPh_list = [] 
    
    Back = Backbone()
    SiPh_list_comp.append(Backbone()) 
    
    #now creating the two benzenes 
    
    #first benzene with updated coordinates 
    benzene1 = coord_add(Benzene(), Back) 
    benzene1 = Translation(benzene1, [0,-1.39,0])
    benzene1_rotated = Rotation(benzene1,90-35.25,-60,0)
    benzene1_rotated = Translation(benzene1_rotated,[0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    SiPh_list_comp.append(benzene1_rotated)
    
    benzene2 = coord_add(Benzene(),benzene1) 
    benzene2 = Translation(benzene2, [0,-1.39,0])
    benzene2_rotated = Rotation(benzene2,35.25-90,-60,0)
    benzene2_rotated = Translation(benzene2_rotated,[0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     -1 *math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    SiPh_list_comp.append(benzene2_rotated)
    
    
    for comp in SiPh_list_comp:
        
        for coords in comp: 
            
            SiPh_list.append(coords) 
                                      
    
    
    return SiPh_list

#### SiEth segment creator ####--------------------------------------------------------

def SiEth():
    
    '''
    creates a eth 
    '''
    
    #initial lists 
    SiEth_list = [] 
    SiEth_comp = [] 
    
    #now creating compenents, moving them around, and adding them to list 
    
    Back = Backbone() 
    SiEth_comp.append(Back) 
    
    Methal1 = coord_add(Methal(), Back) 
    Methal1 = Rotation(Methal1, 90-35.25, 0, 0) 
    Methal1 = Translation(Methal1, [0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    SiEth_comp.append(Methal1)
    
    Ethyl1 = coord_add(Eth(), Methal1) 
    Ethyl1 = Rotation(Ethyl1, 0,-90,0)
    Ethyl1 = Translation(Ethyl1, [0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     -1 *math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    SiEth_comp.append(Ethyl1)
    
    for comp in SiEth_comp: 
        for values in comp: 
            SiEth_list.append(values) 
            
    return SiEth_list
    
#### SiMe segement creator ####-------------------------------------------------------

def SiMe(): 
    
    '''
    adds a methal with a h group onto the backbone 
    '''
    #initial lists 
    SiMe_list = [] 
    SiMe_comp = [] 
    
    #now creating compenents, moving them around, and adding them to list 
    
    Back = Backbone() 
    H1 = [0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_H,
                                     -math.cos(math.radians(35.25))*Bond_Length.Si_H,'H',3] 
    Back.append(H1)
    SiMe_comp.append(Back) 
    
    Methal1 = coord_add(Methal(), Back) 
    Methal1 = Rotation(Methal1, 90-35.25, 0, 0) 
    Methal1 = Translation(Methal1, [0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    SiMe_comp.append(Methal1)
    
    for comp in SiMe_comp: 
        for values in comp: 
            SiMe_list.append(values) 
            
    return SiMe_list
    


#### Segmentor Random ####--------------------------------------------------------------

class Segements() :
    
    
    BenzeneR = Benzene() 
    BenzeneR = Translation(BenzeneR, [0,-1.39,0])
    BenzeneR = Rotation(BenzeneR,90-35.25,-60,0)
    BenzeneR = Translation(BenzeneR,[0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    
    BenzeneL = Benzene() 
    BenzeneL = Translation(BenzeneL, [0,-1.39,0])
    BenzeneL = Rotation(BenzeneL,-90+35.25,-60,0)
    BenzeneL = Translation(BenzeneL,[0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     -1* math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    MethalR = Methal()
    MethalR = Rotation(MethalR, 90-35.25, 0, 0) 
    MethalR = Translation(MethalR, [0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     math.cos(math.radians(35.25))*Bond_Length.Si_C])
    MethalL = Methal()
    MethalL = Rotation(MethalL, -90+35.25, 0, 0) 
    MethalL = Translation(MethalL, [0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     -1* math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    EthylR = Eth() 
    EthylR = Rotation(EthylR, 0,90,0)
    EthylR = Translation(EthylR, [0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    EthylL = Eth()
    EthylL = Rotation(EthylL, 0,-90,0)
    EthylL = Translation(EthylL, [0,
                                     -1* math.sin(math.radians(35.25))*Bond_Length.Si_C,
                                     -1 *math.cos(math.radians(35.25))*Bond_Length.Si_C])
    
    
    
    
    HR= [[0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_H,
                                     math.cos(math.radians(35.25))*Bond_Length.Si_H,'H',1]] 
    HL= [[0,-1* math.sin(math.radians(35.25))*Bond_Length.Si_H,
                                     -math.cos(math.radians(35.25))*Bond_Length.Si_H,'H',1]]
    
    Back = Backbone() 
    
    
    





def SegRand(): 
    '''
    Creates a random segement to hopefully decrease the overall strain
    
    there is 
    '''
    
    #creating the backbone
    
    Back = Segements.Back
    
    
    #creating the two random structures 
    
    r = random.randint(1,7) 
    l = random.randint(1,7) 
    
    rcomp = None
    lcomp = None 
    
    comp_list = [] 
    coord_list = [] 
    
    comp_list.append(Back)
    
    #creating the if statements 
    
    if r > 0 and r <= 3 : 
        rcomp = coord_add(Segements.BenzeneR, Back)
    elif r > 3 and r <= 5 : 
        rcomp = coord_add(Segements.MethalR, Back)
    elif r == 6 :
        rcomp = coord_add(Segements.EthylR, Back)
    else: 
        rcomp = coord_add(Segements.HR, Back)
        
        
    comp_list.append(rcomp)
    
    if l > 0 and l <= 3 : 
        lcomp = coord_add(Segements.BenzeneL, rcomp) 
    elif l > 3 and l <= 5 : 
        lcomp = coord_add(Segements.MethalL, rcomp)  
    elif l == 6 :
        lcomp = coord_add(Segements.EthylL, rcomp)
    else: 
        lcomp = coord_add(Segements.HL, rcomp)    
        
    comp_list.append(lcomp)
        
        
    
    for comp in comp_list: 
        for coord in comp: 
            coord_list.append(coord)
            
    return coord_list

    


#### Segmentor addition function ####---------------------------------------------------
def Segementor():
    '''
    adds segments together to form blocks 
    '''
    
    Seg_comps = [] 
    Seg_list = []
    
    #first SiPh
    
    SiPh_1 = SiPh()
    Seg_comps.append(SiPh_1)
    
    #second 
    
    SiPh_2 = coord_add(SiPh_1,SiPh_1) 
    SiPh_2 = Translation(SiPh_2, [2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0])
    Seg_comps.append(SiPh_2)
    
    #third 
    
    SiPh_3 = coord_add(SiPh_1,SiPh_2) 
    SiPh_3 = Translation(SiPh_3, [2 *2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0])
    
    Seg_comps.append(SiPh_3)
    
    #first SiEth
    
    SiEth_g = SiEth() 
    SiEth_1 = coord_add(SiEth_g,SiPh_3)
    SiEth_1 = Translation(SiEth_1, [3 *2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0])
    Seg_comps.append(SiEth_1)
    
    #second
    
    SiEth_2 = coord_add(SiEth_g,SiEth_1)
    SiEth_2 = Translation(SiEth_2, [4 *2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0])
    Seg_comps.append(SiEth_2)
    
    #first SiMe 
    
    SiMe_g = SiMe() 
    SiMe_1 = coord_add(SiMe_g, SiEth_2)
    SiMe_1 = Translation(SiMe_1, [5 *2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0])
    
    Seg_comps.append(SiMe_1)
    
    #second 
    
    SiMe_g = SiMe() 
    SiMe_2 = coord_add(SiMe_g, SiMe_1)
    SiMe_2 = Translation(SiMe_2, [6 *2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0])
    
    Seg_comps.append(SiMe_2)
    
    
    

    
    for comps in Seg_comps: 
        
        for coords in comps: 
            
            Seg_list.append(coords)
            
    return Seg_list


#### Rand Block Creator ####-------------------------------------------------------


def BlockRand(number): 
    '''
    takes in the rand segments and outputs a block of the length listed 
    
    each 15000 is the number of atoms in a 1000 segment chain 
    '''
    
    comp_list = [] 
    coord_list = [] 
    
    
    
    #first segment
    first_segment = SegRand() 
    comp_list.append(first_segment)
    
    #creating the first methal end group
    end_group1 = Methal() 
    end_group1 = coord_add(end_group1, comp_list[-1])
    end_group1 = Rotation(end_group1, 0, 180, -35.25-90)
    end_group1 = Translation(end_group1, [-Bond_Length.Si_C*math.cos(math.radians(35.25)),
                                          Bond_Length.Si_C*math.sin(math.radians(35.25)),
                                          0])
    comp_list.append(end_group1)
    
    for i in range(number-1) :
        new_seg = coord_add(SegRand(), comp_list[-1])
        new_seg = Translation(new_seg, [(i+1) *2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0])
        comp_list.append(new_seg) 
    
    end_group2= Methal() 
    end_group2 = coord_add(end_group2, comp_list[-1])
    end_group2 = Rotation(end_group2, 0, 180, -35.25+90)
    end_group2 = Translation(end_group2, [Bond_Length.Si_O*math.cos(math.radians(35.25)) * 2 * (number - .5) + 
                                          Bond_Length.C_O*math.cos(math.radians(35.25)),
                                          Bond_Length.Si_O * math.sin(math.radians(35.25))-
                                          Bond_Length.C_O * math.sin(math.radians(35.25)),
                                          0])
    comp_list.append(end_group2)
    
        
    for comps in comp_list: 
        for coord in comps: 
            coord_list.append(coord) 
            
    return coord_list
        
        
        
    
#### Block add ####--------------------------------------------------------------------

def Block(number):
    
    '''
    takes in a number of blocks to be created and outputs a polymer chain with that many blocks 
    uses the block created by a function
    '''
    
    first_segment = Segementor() 
    
    block_comp = [] 
    block_list = []
    
    block_comp.append(first_segment)
    
    for i in range(0,number-1): 
        
        new_segment = coord_add(first_segment, block_comp[i]) 
        new_segment = Translation(new_segment, [(i+1)*7 *2 * math.cos(math.radians(35.25))*Bond_Length.Si_O,
                                  0,
                                  0])
        block_comp.append(new_segment)
        
    for comp in block_comp: 
        for values in comp: 
            block_list.append(values) 
            
    return block_list

def BlockMult(number,block_type): 
    '''
    takes the blocks and expaneds them first in the z direction then in the y direction 
    the number is the number**2 = how many chains there are for the sake of programming time the chains     
    will all be the same but this can be adjusted later 
    
    there are two block types block and block rand 
    '''
    
    comp_listH = [] 
    coord_listH = [] 
    
    #creating the first segment 
    first_segment = block_type 
    comp_listH.append(first_segment)
    
    for i in range(number-1) :
        
        new_seg = coord_add(first_segment, comp_listH[i])
        new_seg = Translation(new_seg, [0,0, (i+1)*10.5]) 
        
        comp_listH.append(new_seg)
        
    for comp in comp_listH : 
        for coord in comp: 
            coord_listH.append(coord)
            
    comp_listV = [] 
    coord_list = [] 
    
    comp_listV.append(coord_listH) 
            
    for i in range(number-1) : 
        
        new_seg = coord_add(coord_listH, comp_listV[i])
        new_seg = Translation(new_seg, [0,(i+1)*6,0])
        
        comp_listV.append(new_seg)
        
    for comp in comp_listV: 
        for coord in comp: 
            coord_list.append(coord) 
            
            
    return coord_list
    
    

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
        
    max_id = max(atom_id_list)
    last_num = str(max_id) + '\n'
            
    with open(filename,'r') as f: 
        lines = f.readlines() 
        
    lines[0] = last_num 
    
    with open(filename, 'w') as f: 
        f.writelines(lines) 
        

def main():
    
    coords =BlockMult(1, BlockRand(100))
    # print(coords)
    
    location = 'G:/Shared drives/1 Kathy Lu Group/Harrison/reaxff/xyz_files/'
    name = 'Rand_w_endgroup6.xyz'

    xyz_appender(coords,location+name)


if __name__ == '__main__':
    main() 
    # pass
