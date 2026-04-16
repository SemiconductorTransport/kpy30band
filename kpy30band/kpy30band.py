#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 18:04:04 2026

@author: badal.mondal
"""

#==============================================================================
from .src import _initiallize_kp_params, _kp_H_30x30_ZB, _AlloyParams
from .basicfns import _SaveData2File
from .utilities import _plot_bandstr
import numpy as np

#==============================================================================
class k_dot_p(_initiallize_kp_params, _kp_H_30x30_ZB):
    def __init__(self, binaries=['Si', 'Ge'], pseudomorphic_strain:bool=False, 
                 substrate:str|float=None, alloy_crystal_structure:str='zb', 
                 use_this_params:dict=None, alloy_type:str=None, 
                 save_file_dir='.', log_info:str='1'):
        self.save_dir_ = save_file_dir
        self.log_output = str(log_info)
        
        _initiallize_kp_params.__init__(self)
        self._initalize_mater_parameters(binaries=binaries,
                pseudomorphic_strain=pseudomorphic_strain, 
                substrate=substrate, alloy_crystal_structure=alloy_crystal_structure, 
                use_this_params=use_this_params, alloy_type=alloy_type)
    
    def kp_30x30(self, compositions:float|np.ndarray|list=None, kpath_list:str|list='L-G', 
                 nkpoints:int=41, return_eigen_val_vec_both:bool=False,
                 save_data:bool=False, save_file:str='eigenval.h5'):
        
        if isinstance(compositions, (int, float)):
            compositions = np.array([compositions], dtype=float)
        else:
            compositions = np.array(compositions, dtype=float)
            
        self._get_alloy_params(compositions)
        
        if self.alloy_crys_type_ == 'zb':
            kpH = _kp_H_30x30_ZB(kpath_list=kpath_list, nkpoints=nkpoints, 
                                 return_eigen_val_vec_both=return_eigen_val_vec_both)
            eigenvalvecs_k_path = {}
            for ii in range(len(compositions)):
                alloy_params_comp = {key: val[ii] for key, val in self.alloy_params_.items()}
                eigenvalvecs_k_path[f'{compositions[ii]:0.5f}'] = kpH._diagonalize_H_30x30(alloy_params_comp)
        else:
            raise ValueError('Hamiltonian for only ZB structures is implemented so far. Contact developer.')
        
        if save_data:
            _SaveData2File._save_data_2_file(data=eigenvalvecs_k_path, save_dir=self.save_dir_, 
                                             file_name=save_file, print_log=self.log_output)
        return eigenvalvecs_k_path
            
class process_data(_SaveData2File):
    def __init__(self, save_file_dir='.', log_info:str=None):
        self.save_dir_ = save_file_dir
        self.log_output = str(log_info)
        
    def save_data_in_file(self, eigenvalvecs_k_path, save_data:bool=False, 
                          save_file:str='eigenval.h5'):
        self._save_data_2_file(data=eigenvalvecs_k_path, save_dir=self.save_dir_, 
                               file_name=save_file, print_log=self.log_output)
        
    def read_data_from_file(self, data_file_name, read_this_compos:list|str=None,
                            read_this_k_paths:list|str=None):
        # Default: read all bands
        return self._read_data_file(f'{self.save_dir_}/{data_file_name}', 
                                    read_this_compos=read_this_compos,
                                    read_this_k_paths=read_this_k_paths)
    
    def process_ek_data(self, data_file_name, read_this_compos:list|str=None,
                        read_this_k_paths:list|str=None):
        return self._process_ek_data(f'{self.save_dir_}/{data_file_name}', 
                                     read_this_compos=read_this_compos,
                                     read_this_k_paths=read_this_k_paths)
    
class plottings(_plot_bandstr):
    """
    Plotting class for mobilitypy.
    """
    def __init__(self, save_figure_dir='.', log_info:str=None):
        """
        Intializing mobilitypy Plotting class.

        Parameters
        ----------
        save_figure_dir : str, optional
            Directory where to save the figure. The default is current directory.

        """
        self.save_figure_directory = save_figure_dir
        self.log_output = str(log_info)
        _plot_bandstr.__init__(self, save_figure_dir=self.save_figure_directory, log_info=self.log_output)
        
    def plot(self, kpts, bands_energy, special_kpts=None, fig=None, ax=None, save_file_name=None, 
              ymin=None, ymax=None, mode:str= 'bandstructure', annotate_pos=(0,0), 
              annotatetextoffset=(0,-20),title_text:str=None, xaxis_label:str=None, 
              yaxis_label:str=r'E (eV)', ls_spkpt='--', lc_spkpt='gray',
              color='k', line_style='-', color_map='viridis', show_legend:bool=False, 
              show_colorbar:bool=False, colorbar_label:str=None, savefig:bool=False,
              vmin=None, vmax=None, show_plot:bool=True, **kwargs_savefig):
        return self._plot(kpts, bands_energy, special_kpts=special_kpts, fig=fig, ax=ax, 
                                  save_file_name=save_file_name, ymin=ymin, ymax=ymax, mode=mode, 
                                  annotate_pos=annotate_pos, annotatetextoffset=annotatetextoffset,
                                  title_text=title_text, xaxis_label=xaxis_label, yaxis_label=yaxis_label, 
                                  ls_spkpt=ls_spkpt, lc_spkpt=lc_spkpt, color=color, line_style=line_style, 
                                  color_map=color_map, show_legend=show_legend, show_colorbar=show_colorbar, 
                                  colorbar_label=colorbar_label, savefig=savefig, vmin=vmin, vmax=vmax, 
                                  show_plot=show_plot, **kwargs_savefig)