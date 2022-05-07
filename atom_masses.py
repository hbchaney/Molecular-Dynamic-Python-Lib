# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 16:11:39 2022

@author: hbcha
"""

'''
H,1.00794
C,12.0107
Si,28.0855
O,15.9994
'''

masses = {'H':1.00794, 'C':12.0107, 'Si':28.0855, 'O':15.9994 }

methane = {'H': 4,'C': 1,'O' :0, 'Si' : 0}  
ethane = {'H': 6,'C': 2,'O' :0, 'Si' : 0}
ethene = {'H': 4,'C': 2,'O' :0, 'Si' : 0}
ethyne = {'H': 2,'C': 2,'O' :0, 'Si' : 0}
dihydrogen = {'H': 2,'C': 0,'O' :0, 'Si' : 0}
c_monoxide = {'H': 0,'C': 1,'O' :1, 'Si' : 0}
c_dioxide = {'H': 0,'C': 1,'O' :2, 'Si' : 0}
water = {'H': 2,'C': 0,'O' :1, 'Si' : 0}
oxygen = {'H': 0,'C': 0,'O' :2, 'Si' : 0}
Benzene = {'H': 6,'C': 6,'O' :0, 'Si' : 0}

mol_list = [methane,ethane,ethene,ethyne,dihydrogen,c_monoxide,c_dioxide,water,oxygen,Benzene]

mass_list = [] 

for mol in mol_list: 
    mass = mol['H']*masses['H'] + mol['C']*masses['C'] + mol['Si']*masses['Si'] + mol['O']*masses['O']
    mass_list.append(mass) 
    
print(mass_list)

    