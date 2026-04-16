#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 14:45:32 2026

@author: badal.mondal
"""
from ._alloy_params import _AlloyParams
import numpy as np

## ==============================================================================
class _initiallize_kp_params(_AlloyParams):
    def __init__(self, print_log=None):
        self.print_info = print_log
    
    def _initalize_mater_parameters(self, binaries=['Si', 'Ge'], pseudomorphic_strain:bool=False, 
                                    substrate:str|float=None, alloy_crystal_structure:str='zb', 
                                    use_this_params:dict=None, alloy_type:str=None):
        if pseudomorphic_strain and (substrate is None):
            # This allows to return Error in the very begining without starting any calculations.
            raise ValueError('substrate tag can not be None when pseudomorphic_strain=True.')

        _AlloyParams.__init__(self, binaries=binaries, 
                              alloy_crystal_structure=alloy_crystal_structure,
                              alloy_type=alloy_type, use_this_params=use_this_params)
        
        if pseudomorphic_strain: 
            self.strain_tensor = self._cal_pseudomorphic_strain(substrate)
        
#%% ===========================================================================
class _kpoints_from_path:
    kpath_map = {'G':[0,0,0], 'L':[0.5,0.5,0.5], 'X':[0,0.5,0.5]}
    def __init__(self):
        pass
    @classmethod
    def _generate_k_points_from_path(cls, kpath='L-G', nkpts=41):
        # In reciprocal space
        kpath = [kpath] if isinstance(kpath, str) else kpath
        k_points_dict = {}
        for kpath_frac in kpath:
            k_point_name = kpath_frac.split('-')
            assert len(k_point_name)==2, 'Single section is allowed. E.g. G-L'
            k_points_dict[kpath_frac] = np.linspace(cls.kpath_map[k_point_name[0]], 
                                                    cls.kpath_map[k_point_name[1]], 
                                                    nkpts)
        return k_points_dict  