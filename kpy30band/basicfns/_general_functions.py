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
            # Data structure: composition group => k-path sheets
            with h5py.File(fname_save_file, 'w') as f:
                for compos, k_path_vals in data.items():
                    data_group = f.create_group(compos)
                    for k_path, vals in k_path_vals.items():
                        data_group.create_dataset(k_path, data=vals, dtype='d')
        return fname_save_file
    
    @classmethod
    def _save_data_2_file(cls, data=None, save_dir='.', file_name:str='', 
                          print_log:bool=False, print_msg:str='Saving to file...'):
        header_txt = '! kx(nm^-1)  ky(nm^-1)  kz(nm^-1) E(eV)'          
            
        if print_log > '1': print(f"{'='*cls._draw_line_length}\n- {print_msg}.")
        
        fname_save_file = cls._save_2_file(data=data, save_dir=save_dir, 
                                           file_name=file_name, 
                                           header_txt=header_txt, 
                                           footer_txt='',
                                           comments_symbol='!')
        if print_log > '1': print(f'-- Filepath: {fname_save_file}\n- Done')
        return fname_save_file
    
        
    @staticmethod
    def _read_data_file(full_file_path_name:str='', read_this_compos:list|str=None,
                        read_this_k_paths:list|str=None):
        plot_data = {}
        with h5py.File(full_file_path_name, 'r') as f:
            # Data structure: composition group => k-path sheets
            read_compos = list(f.keys()) if read_this_compos is None else read_this_compos
            if read_this_k_paths is None:
                read_f = f[read_compos[0]].keys()
            elif isinstance(read_this_k_paths, str):
                read_f = [read_this_k_paths]
            else:
                read_f = read_this_k_paths
            for compos in read_compos: # composition loop
                plot_data[compos] = {}
                for k_paths in read_f: # k-path loop
                    plot_data[compos][k_paths] = f[compos][k_paths][:]
        return plot_data
    
    @classmethod
    def _process_ek_data(cls, data_file_name, read_this_compos:list|str=None,
                         read_this_k_paths:list|str=None):
        # Default: read all bands
        read_data = cls._read_data_file(data_file_name,
                                        read_this_compos=read_this_compos,
                                        read_this_k_paths=read_this_k_paths) 
        for_compos = {}
        for compos in read_data.keys(): # loop over compositions
            shift = 0
            special_kpts_pos, special_kpts_label, kpts, bands_e = [], [], [], []
            for k_paths, r_data in read_data[compos].items(): # loop over k-path
                k_components = k_paths.split('-')
                k_ps = r_data[:, 0:3] # First 3 numbers are k-points
                XX = np.sqrt(np.sum(k_ps*k_ps, axis=1))
                #-----------------------------------------------------
                # This is to move x-axis coordinates to make the plot 
                if len(k_components) > 1:
                    XX_fi = XX[-1] - XX[0]
                    XX = XX + shift if XX_fi > 0 else XX[0]-XX+shift
                    shift = XX[-1]
                    special_kpts_pos += [XX[0], XX[-1]]
                else:
                    special_kpts_pos += [XX[0]]
                #-----------------------------------------------------
                kpts.append(XX)
                bands_e.append(r_data[:, 3:].T) # 3 is for kx, ky, kz at the begining
                special_kpts_label += k_components   
            kpts = np.concatenate(kpts, axis=0)
            bands_e = np.concatenate(bands_e, axis=1)
            for_compos[compos] = {'shifted_k(nm^-1)':kpts, 'E(ev)': bands_e, 
                                  'special_kpts': (special_kpts_pos, special_kpts_label)}
        return for_compos
    