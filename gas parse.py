# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 21:29:01 2021

@author: Harrison
"""

deleted_atoms = '420 421 453 7349 422 743 774 804 944 955 956 976 2909 1077 1078 1263 1265 1429 2223 2254 2284 2424 2312 11401 2435 2436 2521 4170 2557 2558 2743 2745 2837 8829 3380 3381 3413 3382 3703 3734 3764 3904 3861 6961 3915 3916 3923 10309 4037 4038 4223 4225 4389 8896 4860 4861 5717 4893 4862 5183 5214 5244 5384 5395 5396 5481 9530 5486 5488 5487 5502 5612 5517 5518 5703 5705 5869 6663 6694 6724 6864 6875 6876 6966 6968 6982 6967 7092 7117 6997 6998 7035 7037 7115 7036 7124 7183 7185 8143 8174 8204 8344 8355 8356 8477 8478 8663 8665 9623 9654 9684 9824 9835 9836 9921 11493 9957 9958 10143 10145 11103 11134 11164 11304 11315 11316 11437 11438 11623 11625'

deleted_atoms = deleted_atoms.split(' ')
# print(deleted_atoms)

deleted_atoms = [int(x) for x in deleted_atoms]
# print(deleted_atoms)

with open('rand_2x2_t4_500000.in') as f: 
    lines = f.readlines() 
    
# print(lines[16:])
lines = lines[16:]

# print(lines[-1])
atom_coord = {}
for stuff in lines:
    stuff = stuff.split(' ')
    atom_id = int(stuff[0])
    atom_type = int(stuff[1])
    
    atom_coord[atom_id] = atom_type

# print(atom_coord)

for atoms in deleted_atoms:
    print(atom_coord[atoms])
    