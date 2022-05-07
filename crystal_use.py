# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 09:38:41 2022

@author: hbcha
"""
#################################################################################
#######################################TESTS#####################################
#################################################################################
from Crystal_parser import bond_frame,loc_frame,combined_frame
# x.file_read(r'G:\My Drive\Research\reaxff\VMHM\VMHM_t2_gas\VMHM_t2_9\gas_creation.51.reaxc')


# x = combined_frame() 
# x.file_read(r'G:\My Drive\Research\reaxff\VMHM\VMHM_t2_gas\VMHM_t2_9\gas_creation.51.reaxc','')
# print(x.bond_data)

# x = loc_frame() 
# x.file_read(r'G:\My Drive\Research\reaxff\VMHM\VMHM_t2_9\loop_heat.51')

x = combined_frame() 
x.file_read(r'G:\My Drive\Research\reaxff\VMHM\VMHM_t2_gas\VMHM_t2_9\gas_creation.51.reaxc', r'G:\My Drive\Research\reaxff\VMHM\VMHM_t2_9\loop_heat.51')
x.comb_del(x.bond_exclude(1,3))
x.comb_del(x.bond_exclude(1,1))
x.comb_del(x.bond_exclude(1,2,reverse=1))


print(x)
x.to_file(r'G:\My Drive\Research\reaxff\VMHM\VMHM_crystal_info\Si_O_noSi_Si.in')


