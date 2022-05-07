# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 11:54:13 2022

@author: hbcha
"""
def wt(Si,O,C,H):
    
    Si_m = 28.085
    O_m = 15.999
    C_m = 12.011
    H_m = 1.008
    
    at_t = Si + O + C + H 
    
    wt_t = Si*Si_m+O*O_m+C*C_m+H*H_m
    
    Siat, Oat, Cat, Hat = Si/at_t,O/at_t,C/at_t,H/at_t
    Siwt,Owt,Cwt,Hwt = Si*Si_m/wt_t,O*O_m/wt_t,C*C_m/wt_t,H*H_m/wt_t
    
    print('atom %-')
    print(f'Si : {round(Siat*100,2)}')
    print(f'O : {round(Oat*100,2)}')
    print(f'C : {round(Cat*100,2)}')
    print(f'H : {round(Hat*100,2)}')
    
    print('weight %-')
    print(f'Si : {round(Siwt*100,2)}')
    print(f'O : {round(Owt*100,2)}')
    print(f'C : {round(Cwt*100,2)}')
    print(f'H : {round(Hwt*100,2)}')
    
    
wt(1,1,2,6)

    