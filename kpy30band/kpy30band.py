#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 18:04:04 2026

@author: badal.mondal
"""

#==============================================================================
from .src import _kp_H_30x30

#==============================================================================
class k_dot_p:
    def __init__(self):
        pass
    
    def kp_30x30(self, a0:float, E_Gamma:dict, E_Delta:dict, Off_diag:dict, 
                 kpath_list:str|list='L-G', nkpoints:int=41, 
                 return_eigen_val_vec_both:bool=False):
        _kp_H_30x30.__init__(a0, E_Gamma, E_Delta, Off_diag, kpath_list=kpath_list, nkpoints=nkpoints, 
                             return_eigen_val_vec_both=return_eigen_val_vec_both)
        return self._diagonalize_H_30x30()
        