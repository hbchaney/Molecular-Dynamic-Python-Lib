# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 10:14:11 2022

@author: hbcha
"""


import Harrison 

# x1 = Harrison.linear_Polymer.Si_ORandCreator(112, {'Phenyl':0,'Methyl':1,'Vinyl':1,'Hydride':0})
# x2 = Harrison.linear_Polymer.Si_ORandCreator(112, {'Phenyl':0,'Methyl':1,'Vinyl':1,'Hydride':0})
# x3 = Harrison.linear_Polymer.Si_ORandCreator(112, {'Phenyl':0,'Methyl':1,'Vinyl':1,'Hydride':0})
# x4 = Harrison.linear_Polymer.Si_ORandCreator(112, {'Phenyl':0,'Methyl':1,'Vinyl':1,'Hydride':0})
# x5 = Harrison.linear_Polymer.Si_ORandCreator(112, {'Phenyl':0,'Methyl':1,'Vinyl':1,'Hydride':0})
# x6 = Harrison.linear_Polymer.Si_ORandCreator(112, {'Phenyl':0,'Methyl':1,'Vinyl':1,'Hydride':0})
# x7 = Harrison.linear_Polymer.Si_ORandCreator(112, {'Phenyl':0,'Methyl':1,'Vinyl':1,'Hydride':0})
# x8 = Harrison.linear_Polymer.Si_ORandCreator(28, {'Phenyl':0,'Methyl':1,'Vinyl':1,'Hydride':0})

# y1 = Harrison.linear_Polymer.Si_ORandCreator(108, {'Phenyl':0,'Methyl':1,'Vinyl':0,'Hydride':1})
# y2 = Harrison.linear_Polymer.Si_ORandCreator(108, {'Phenyl':0,'Methyl':1,'Vinyl':0,'Hydride':1}) 

# # print(x1,x2,x3)


# z = x1 + x2 + x3 + x4 + y1 + y2 + x5 + x6 + x7 + x8 

# z.to_file('SAX_VMHM_wt.in',r'G:\My Drive\Research\reaxff\lammps_in\Blends')


x1 = Harrison.linear_Polymer.Si_ORandCreator(133, {'Phenyl':1,'Methyl':1,'Vinyl':0,'Hydride':0})
x2 = Harrison.linear_Polymer.Si_ORandCreator(133, {'Phenyl':1,'Methyl':1,'Vinyl':0,'Hydride':0})
x3 = Harrison.linear_Polymer.Si_ORandCreator(133, {'Phenyl':1,'Methyl':1,'Vinyl':0,'Hydride':0})
x4 = Harrison.linear_Polymer.Si_ORandCreator(133, {'Phenyl':1,'Methyl':1,'Vinyl':0,'Hydride':0})
x5 = Harrison.linear_Polymer.Si_ORandCreator(133, {'Phenyl':1,'Methyl':1,'Vinyl':0,'Hydride':0})
x6 = Harrison.linear_Polymer.Si_ORandCreator(133, {'Phenyl':1,'Methyl':1,'Vinyl':0,'Hydride':0})

y1 = Harrison.linear_Polymer.Si_ORandCreator(108, {'Phenyl':0,'Methyl':1,'Vinyl':0,'Hydride':1})
y2 = Harrison.linear_Polymer.Si_ORandCreator(108, {'Phenyl':0,'Methyl':1,'Vinyl':0,'Hydride':1})
y3 = Harrison.linear_Polymer.Si_ORandCreator(108, {'Phenyl':0,'Methyl':1,'Vinyl':0,'Hydride':1}) 



z = x1 + x2 + x3  + y1 + y2 + x4 + x5 + x6 +y3 

print(z)

z.to_file('SAX_MPHM_wt_t2_1.in',r'G:\My Drive\Research\reaxff\lammps_in\Blends')