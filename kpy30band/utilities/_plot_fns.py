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
    def __init__(self, save_figure_dir='.', log_info=None):
        """
        Initialize the plotting class.

        Parameters
        ----------
        save_figure_dir : str/path, optional
            Directory where to save the figure. The default is current directory.
       """
        _GeneratePlots.__init__(self, save_figure_dir=save_figure_dir, log_info=log_info)
        
    def _plot(self, kpts, bands_energy, special_kpts=None, fig=None,  
              ax=None, save_file_name=None, ymin=None, ymax=None, 
              annotate_text={'text':None, 'pos':(0,0)}, 
              y_axis_major_tick:float=None, title_text:str=None,  
              xaxis_label:str=None, yaxis_label:str=r'E (eV)', ls_spkpt='--', 
              lc_spkpt='gray', line_marker='o', line_style='-', 
              color='gray', color_map='viridis', show_legend:bool=False, 
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

        if title_text is not None: self.ax.set_title(title_text)
        if annotate_text.get('text'): 
            self.ax.annotate(annotate_text['text'], (0,0), 
                             xytext=annotate_text['pos'], 
                             textcoords='axes fraction',
                             va="center", ha="center")
        
        #ax, return_plot = self._plot_2d_plane(results, ax, color=color, ls=ls_2d)
        if self.plot_mode == 'band_structure':
            return_plot = self._plot_bandstructure(kpts, bands_energy, color=color, 
                                                   ls=line_style,marker=line_marker)
        else:
            raise ValueError(f'Requested plot mode {self.plot_mode} has not been implemented yet. Contact developer.')
            
        if show_colorbar and (self.fig is not None):
            cbar = self.fig.colorbar(return_plot)
            if colorbar_label is not None:
                cbar.set_label(colorbar_label)
                
        self.ax.axhline(y=0, color=lc_spkpt, ls=ls_spkpt)
        if special_kpts is not None:
            xticks_ = _plot_bandstr._map_special_kpoints_labels(special_kpts[0], special_kpts[1])
            xticks_pos = np.array(list(xticks_.keys()), dtype='float')
            self.ax.set_xticks(xticks_pos, labels=list(xticks_.values()))
            for xx in xticks_pos:
                self.ax.axvline(x=xx, color=lc_spkpt, ls=ls_spkpt)
        if y_axis_major_tick: self.ax.yaxis.set_major_locator(ticker.MultipleLocator(y_axis_major_tick))
        if save_file_name is None:
            if show_plot: plt.show()
        else:
            self._save_figure(save_file_name, savefig=savefig, show_plot=show_plot, 
                              fig=self.fig, **kwargs_savefig)
            
        return 
    
    @classmethod
    def _map_special_kpoints_labels(cls, special_kpts_pos, special_kpts_label):
        _name_map = {'G': r'$\Gamma$', 'D': r'$\Delta$', 'Z_romb': 'Z', 
                     'F_romb': 'F', 'L_romb': 'L', 'Z_tet': 'Z', 'X_tet': 'X',
                     'N_tet':'N'}
        xticks_ = {}
        for ii, pos_ in enumerate(special_kpts_pos):
            key = f'{pos_:0.3f}'
            ll = _name_map.get(special_kpts_label[ii])
            if ll is None: ll = special_kpts_label[ii]
            if key in xticks_:
                if xticks_[key] == ll: continue
                xticks_[key] += f',{ll}'
            else:
                xticks_[key] = ll
        return xticks_
        
    def _plot_bandstructure(self, kpts, bands_e, color=None, ls='-', marker='o'):
        for be in bands_e:
            self.ax.plot(kpts, be, ls=ls, color=color, marker=marker)
        return
        
    