#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 18:04:04 2026

@author: badal.mondal
"""

#==============================================================================
from .src import _DataBase, _initiallize_kp_params, _kp_H_30x30_ZB, _AlloyParams
from .basicfns import _SaveData2File
from .utilities import _plot_bandstr
import numpy as np

## ==============================================================================
class DataBase(_DataBase):
    """
    The functions in this class use to print/manipulate the material parameters
    in the database.
    """
    def __init__(self):
        self.dtbase = _DataBase()
        
    def print_database(self, for_material:str=None):
        """
        This function prints the information of material parameters from the database.

        Parameters
        ----------
        for_material : string (case sensitive), optional
            The material name for which the parameters will be printed. 
            The name should match the name in database. If None, prints general
            information about the database.
            The default is None.

        Returns
        -------
        None.

        """
        self.dtbase._print_database(for_material=for_material)
    
    def print_materials_available_in_database(self):
        """
        Print the materials available so far in the database.

        Returns
        -------
        None.

        """
        self.dtbase._print_materials_available_in_database()

#==============================================================================
class k_dot_p(_initiallize_kp_params, _kp_H_30x30_ZB):
    def __init__(self, binaries=['Si', 'Ge'], pseudomorphic_strain:bool=False, 
                 substrate:str|float=None, growth_direction_hkl:list[int]=[0,0,1], 
                 alloy_crystal_structure:str='zb', 
                 use_this_params:dict=None, alloy_type:str=None, 
                 save_file_dir='.', log_info:str='1'):
        self.save_dir_ = save_file_dir
        self.log_output = str(log_info)
        
        _initiallize_kp_params.__init__(self)
        self._initalize_mater_parameters(binaries=binaries,
                        pseudomorphic_strain=pseudomorphic_strain, 
                        substrate=substrate, growth_direction=growth_direction_hkl, 
                        alloy_crystal_structure=alloy_crystal_structure, 
                        use_this_params=use_this_params, alloy_type=alloy_type)
    
    def kp_30x30(self, compositions:float|np.ndarray|list=None,
                 overwrite_strain:float|list|np.ndarray=None,
                 kpoints_list:str|list='L-G', nkpoints:int=41, 
                 return_eigen_val_vec_both:bool=False,
                 save_data:bool=False, save_file:str='eigenval.h5'):
        
        # overwrite_strain is the strain values for each composition. 
        # If len(overwrite_strain) > len(compositions): only first len(compositions) values
        # will be used for overwrite_strain, i.e. overwrite_strain_new = overwrite_strain[:len(compositions)]
        # Note: This does not run array of strain values for a single composition number.
        # In that case, one needs to explicitely loop over strain values outside this package.
        
        if isinstance(compositions, (int, float)):
            compositions = np.array([compositions], dtype=float)
        else:
            compositions = np.array(compositions, dtype=float)
            
        self._get_alloy_params(compositions)
        
        strain_tensor = (self._cal_pseudomorphic_strain(overwrite_strain)).T \
                                if self.apply_strain_ else [None]*len(compositions)

        if self.alloy_crys_type_ == 'zb':
            _kp_H_30x30_ZB.__init__(self, kpoints_list=kpoints_list, nkpoints=nkpoints, 
                                    return_eigen_val_vec_both=return_eigen_val_vec_both)
            eigenvalvecs_k_path = {}
            for ii in range(len(compositions)):
                self.strain_tensor = strain_tensor[ii] 
                #print(self.strain_tensor)
                alloy_params_comp = {key: val[ii] for key, val in self.alloy_params_.items()}
                eigenvalvecs_k_path[f'{compositions[ii]:0.5f}'] = self._diagonalize_H_30x30(alloy_params_comp)
        else:
            raise ValueError('Hamiltonian for only ZB structures is implemented so far. Contact developer.')
        
        if save_data:
            _SaveData2File._save_data_2_file(data=eigenvalvecs_k_path, save_dir=self.save_dir_, 
                                             file_name=save_file, print_log=self.log_output)
        return eigenvalvecs_k_path

#==============================================================================           
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

#==============================================================================
class plottings(_plot_bandstr):
    """
    Plotting class.
    """
    def __init__(self, save_figure_dir='.', log_info:str=None):
        """
        Intializing Plotting class.

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
#==============================================================================