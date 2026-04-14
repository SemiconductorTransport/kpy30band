#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 08:42:02 2026

@author: badal.mondal
"""
import numpy as np
import h5py
## ============================================================================
class _SaveData2File:
    _draw_line_length = 72
    def __init__(self):
        pass
    
    @staticmethod
    def _save_2_file(data=None, save_dir='.', file_name:str='',  
                     header_txt:str='', footer_txt:str='',comments_symbol='! ',
                     np_data_fmt='%12.8f', print_log:bool=True):
        if data is None: return
        if not file_name.endswith(('.h5','.dat')):
            raise ValueError('Only .h5 and .dat files are allowed.')
        
        fname_save_file = f'{save_dir}/{file_name}'
        if file_name.endswith('.dat'):
            raise ValueError('Only .h5 files are allowed at the moment.')
            with open(fname_save_file, 'w') as f:
                np.savetxt(f, data, header=header_txt, footer=footer_txt, 
                           fmt=np_data_fmt, comments=comments_symbol)
        else:
            with h5py.File(fname_save_file, 'w') as f:
                for k_path, vals in data.items():
                    f.create_dataset(k_path, data=vals, dtype='d')
        return fname_save_file
    
    @classmethod
    def _save_data_2_file(cls, data=None, save_dir='.', file_name:str='', 
                          print_log:bool=False, print_msg:str='Saving to file...'):
        header_txt = '! kx(2pi/a)  ky(2pi/a)  kz(2pi/a) E(eV)'          
            
        if print_log > '1': print(f"{'='*cls._draw_line_length}\n- {print_msg}.")
        
        fname_save_file = cls._save_2_file(data=data, save_dir=save_dir, 
                                           file_name=file_name, 
                                           header_txt=header_txt, 
                                           footer_txt='',
                                           comments_symbol='!')
        if print_log > '1': print(f'-- Filepath: {fname_save_file}\n- Done')
        return fname_save_file
    
        
    @staticmethod
    def _read_data_file(full_file_path_name:str='', read_this_k_paths:list|str=None):
        plot_data = {}
        with h5py.File(full_file_path_name, 'r') as f:
            read_f = f.keys() if read_this_k_paths is None else read_this_k_paths
            for k_paths in read_f:
                plot_data[k_paths] = f[k_paths][:]
        return plot_data