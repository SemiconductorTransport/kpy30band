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
                                    substrate:str|float='test_sample', growth_direction:list[int]=[0,0,1],
                                    alloy_crystal_structure:str='zb', 
                                    use_this_params:dict=None, alloy_type:str=None):
        if pseudomorphic_strain and (substrate is None):
            # This allows to return Error in the very begining without starting any calculations.
            raise ValueError('Substrate can not be None when pseudomorphic_strain=True.')

        _AlloyParams.__init__(self, binaries=binaries, 
                              alloy_crystal_structure=alloy_crystal_structure,
                              alloy_type=alloy_type, use_this_params=use_this_params)
        
        if pseudomorphic_strain: 
            self.apply_strain_ = True
            self.biaxial_substrate = substrate
            self.growth_hkl = np.array(growth_direction, dtype=float)
            
        
#%% ===========================================================================
class _kpoints_from_point_path:
    def __init__(self, crys_type_):
        if crys_type_ in ['zb', 'dm', 'cube']:
            self.kpath_map = {'G':[0.0,0.0,0.0], 'L':[0.5,0.5,0.5], 'X':[0.0,0.0,1.0], 
                              'K':[0.75,0.75,0.0], 'W':[0.5,0.0,1.0], 'U':[0.25,0.25,1.0]}
        else:
            raise ValueError(f'Special k-points for {crys_type_} is not available yet. Contact developer.')

    def _generate_k_points_from_point_path(self, kpath:str|list='L-G', nkpts:int=11):
        ## kpath = 'L-G'; ['L-G']; ['L-G','G-X']; [L, G, X]; [0,0.5,0.5];
        ##         [[0,0.5,0.5], [0,0.45,0.45], [0,0.35,0.35], ...];
        ##         [[0,0.5,0.5], [0,0.45,0.45], 'L-G', 'L', ...] # mixed array
        # In reciprocal space
        k_points_dict = {}
        
        if isinstance(kpath, str):
            k_points_dict[kpath] = self._generate_k_points_from_k_str(kpath, nkpts)
        else:
            try: # pure numbers or pure strings
                kpath = np.array(kpath)
                if kpath.ndim == 1 and all(isinstance(x, (int, float)) for x in kpath):
                    return {'1': np.array([kpath])}
                elif kpath.ndim == 2:
                    return {'1': np.array(kpath, dtype=float)} # 2d array of floats
                else:
                    pass
            except: # mixed array
                pass
            for ii, kpath_frac in enumerate(kpath):
                if isinstance(kpath_frac, str):
                    k_points_dict[kpath_frac] = self._generate_k_points_from_k_str(kpath_frac, nkpts)
                else:
                    k_points_dict[str(ii)] =  np.array([kpath_frac])
        return k_points_dict  
    
    def _generate_k_points_from_k_str(self, kpath_frac, nkpts):
        k_point_name = kpath_frac.split('-') 
        lenkpts = len(k_point_name)
        if lenkpts == 2:
            for kname in k_point_name:
                self._check_special_kp_name_in_defined(kname)
            return np.linspace(self.kpath_map[k_point_name[0]], 
                               self.kpath_map[k_point_name[1]], nkpts)
        elif lenkpts > 2:
            raise ValueError('only two kpoint name is allowed in k_point_name. E.g., L-G')
        else:
            self._check_special_kp_name_in_defined(k_point_name[0])
            return np.array([self.kpath_map[k_point_name[0]]])
        
    def _check_special_kp_name_in_defined(self, kname):
        if kname not in self.kpath_map.keys():
            raise ValueError(f'{kname} special k-point is not defined in database yet. contact developer.')