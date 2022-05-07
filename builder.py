# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 11:09:22 2022

@author: hbcha
"""

import Harrison 



# x = Harrison.linear_Polymer.Si_ORandCreator(25,{'Phenyl':0,'Methyl':1,'Vinyl':0,'Hydride':1})


# x = Harrison.amorphous_box.file_read(r'G:\My Drive\Research\reaxff\PDMS\PDMS_t2_3','loop_heat.1',[1,2,3,4])
# x = x.rand_expansion((0,0,0),(1,0,0),buffer=0)
# x = x.rand_expansion((0,0,0),(0,1,0),buffer=0)
# x = x.rand_expansion((0,0,0),(0,0,1),dup_num=2,buffer=0)
# x.to_file('PDMS_t4_1.in','G:/My Drive/Research/reaxff/lammps_in/PDMS')

# x = Harrison.amorphous_box.file_read(r'G:\My Drive\Research\reaxff\PHMS\lc_stab_t1_4','loop_heat.38',[4,3,2,1])
# x = x.rand_expansion((0,0,0),(1,0,0),dup_num=2)
# x = x.rand_expansion((0,0,0),(0,1,0))
# x = x.rand_expansion((0,0,0),(0,0,1),dup_num=2)
# x.to_file('PHMS_expanded_1.in','G:/My Drive/Research/reaxff/lammps_in/PHMS')

# print(x)    

# x = Harrison.amorphous_box.file_read(r'G:\My Drive\Research\reaxff\VMHM\VMHM_t2_3','loop_heat.1',[1,2,3,4])
# x = x.rand_expansion((0,0,0),(1,0,0),buffer=0)
# x = x.rand_expansion((0,0,0),(0,1,0),buffer=0)
# x = x.rand_expansion((0,0,0),(0,0,1),buffer=0)
# x.to_file('VMHM_t3_1.in','G:/My Drive/Research/reaxff/lammps_in/VMHM')

# x = Harrison.amorphous_box.file_read(r'G:\My Drive\Research\reaxff\SAX\MPHM_t2_3','loop_heat.1',[1,2,3,4])
# x = x.rand_expansion((0,0,0),(1,0,0),buffer=0,dup_num=2)
# x = x.rand_expansion((0,0,0),(0,1,0),buffer=0)
# # x = x.rand_expansion((0,0,0),(0,0,1),buffer=0)
# x.to_file('MPHM_t3_1.in','G:/My Drive/Research/reaxff/lammps_in/MPHM')

# x = Harrison.amorphous_box.file_read(r'G:\My Drive\Research\reaxff\PHMS\lc_stab_t1_4',r'PDES_t1_3.rdx')
# x = x.rand_expansion((0,0,0),(1,0,0),buffer=0)
# x = x.rand_expansion((0,0,0),(0,1,0),buffer=0)
# x = x.rand_expansion((0,0,0),(0,0,1),buffer=0)
# x.to_file('PDES_t1_4.in',r'G:\My Drive\Research\reaxff\lammps_in\PDES')

# x = Harrison.amorphous_box.file_read(r'G:\My Drive\Research\reaxff\SQHM','SQHM_t1_1.rdx',[1,2,3,4])
# x = x.rand_expansion((0,0,0),(1,0,0),buffer=0)
# x = x.rand_expansion((0,0,0),(0,1,0),buffer=0)
# x = x.rand_expansion((0,0,0),(0,0,1),buffer=0)

x = Harrison.amorphous_box.file_read(r'G:\My Drive\Research\reaxff\PDES\PDES_t1_4', 'loop_heat.4')
x.to_file('PDES_t2_1.in',r'G:\My Drive\Research\reaxff\lammps_in\PDES')


