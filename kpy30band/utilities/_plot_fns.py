#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 14:43:37 2026

@author: badal.mondal
"""

import numpy as np
import matplotlib.pyplot as plt
from ._general_plot_fns import _GeneratePlots
import matplotlib.ticker as ticker

### ===========================================================================
class _plot_bandstr(_GeneratePlots):
    """
    The functions in this class plots the different mobility figures.

    """
    _special_kp_map = {'G':r'$\Gamma$', 'L':'L', 'X':'X', 'G':r'$\Gamma$'}
    def __init__(self, save_figure_dir='.', log_info=None):
        """
        Initialize the plotting class.

        Parameters
        ----------
        save_figure_dir : str/path, optional
            Directory where to save the figure. The default is current directory.
       """
        _GeneratePlots.__init__(self, save_figure_dir=save_figure_dir, log_info=log_info)
        
    def _plot(self, kpts, bands_energy, special_kpts=None, fig=None, ax=None, save_file_name=None, 
              ymin=None, ymax=None, mode:str= 'bandstructure', annotate_pos=(0,0), 
              annotatetextoffset=(0,-20),title_text:str=None, xaxis_label:str=None, 
              yaxis_label:str=r'E (eV)', ls_spkpt='--', lc_spkpt='gray',
              color='gray', line_style='-', color_map='viridis', show_legend:bool=False, 
              show_colorbar:bool=False, colorbar_label:str=None, savefig:bool=False,
              vmin=None, vmax=None, show_plot:bool=True, **kwargs_savefig):
        if ax is None: 
            self.fig, self.ax = plt.subplots(constrained_layout=True)
        else:
            self.fig = fig; self.ax = ax
            
        if yaxis_label is None: yaxis_label=''
        if xaxis_label is None: xaxis_label=''

        self.ax.set_ylabel(yaxis_label)
        self.ax.set_xlabel(xaxis_label)
        self.ax.set_ylim([ymin, ymax])
        self.ax.set_xlim([min(kpts), max(kpts)])

        if title_text is not None: ax.set_title(title_text)
        
        #ax, return_plot = self._plot_2d_plane(results, ax, color=color, ls=ls_2d)
        if mode == 'bandstructure':
            self._plot_bandstructure(kpts, bands_energy, special_xticks=special_kpts, 
                                     color=color, ls=line_style, ls_spkpt=ls_spkpt, 
                                     lc_spkpt=lc_spkpt)
            return_plot = None
        else:
            raise ValueError('Requested plot mode is not allowed')
            
        if show_colorbar and (self.fig is not None):
            cbar = self.fig.colorbar(return_plot)
            if colorbar_label is not None:
                cbar.set_label(colorbar_label)
                
        self.ax.axhline(y=0, color=lc_spkpt, ls=ls_spkpt)
        if special_kpts is not None:
            xticks_pos = special_kpts.values()
            xticks_label = [self._special_kp_map[kk] for kk in special_kpts.keys()]  
            self.ax.set_xticks(list(xticks_pos), labels=list(xticks_label))
            for xx in xticks_pos:
                self.ax.axvline(x=xx, color=lc_spkpt, ls=ls_spkpt)

        if save_file_name is None:
            if show_plot: plt.show()
        else:
            self._save_figure(save_file_name, savefig=savefig, show_plot=show_plot, 
                              fig=self.fig, **kwargs_savefig)
            
        return self.fig, self.ax
    
    def _plot_bandstructure(self, kpts, bands_e, special_xticks=None, color=None, ls='-',
                            ls_spkpt='--', lc_spkpt='gray'):
        for be in bands_e:
            self.ax.plot(kpts, be, ls=ls, color=color)
        return
        
    