#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 14:39:21 2026

@author: badal.mondal
"""

import matplotlib.pyplot as plt
from pathlib import Path
### ===========================================================================

class _GeneratePlots:
    def __init__(self, save_figure_dir='.', log_info=None):
        """
        Initialize the plotting class.

        Parameters
        ----------
        save_figure_dir : str/path, optional
            Directory where to save the figure. The default is current directory.

        Returns
        -------
        None.

        """
        self.log_ = log_info
        self.save_figure_folder = save_figure_dir
        params = {'figure.figsize': (8, 6),
                  'legend.fontsize': 18,
                  'legend.title_fontsize': 18,
                  'axes.labelsize': 24,
                  'axes.titlesize': 24,
                  'xtick.labelsize':24,
                  'xtick.major.width':2,
                  'xtick.major.size':5,
                  'xtick.minor.width':2,
                  'xtick.minor.size':3,
                  'ytick.labelsize': 24,
                  'ytick.major.width':2,
                  'ytick.major.size':5,
                  'ytick.minor.width':2,
                  'ytick.minor.size':3,
                  'errorbar.capsize':2,
                  'lines.markersize':12,
                  'lines.linewidth':2, 
                  'lines.linestyle':'-'}
        plt.rcParams.update(params)
        plt.rc('font', size=24)

    def _save_figure(self,fig_name, savefig:bool=True, show_plot:bool=True,
                     fig=None,  **kwargs_savefig):
        """
        The function saves the matplotlib figure.

        Parameters
        ----------
        fig_name : str
            name of the figure to use including figure extension.
        savefig : bool, optional
            Save the figure or not. The default is True.
        show_plot : bool, optional
            Show the plot or not. The default is True.
        fig : matplotlib figure instance, optional
            matplotlib figure instance. The default is None.
        **kwargs_savefig : matplotlib savefig kwargs
            Any other matplotlib savefig keywords arguments.

        Returns
        -------

        """
        if not savefig: 
            if show_plot: plt.show()
            return 

        Path(self.save_figure_folder).mkdir(parents=True, exist_ok=True)
        save_path = f'{self.save_figure_folder}/{fig_name}'
        if self.log_: print(f'* Plot saved: {save_path}')
        if fig is not None:
            fig.savefig(save_path, bbox_inches='tight', **kwargs_savefig)
        else:
            plt.savefig(save_path, bbox_inches='tight', **kwargs_savefig)
        plt.close()
        return 

    def save_figure(self,fig_name, savefig:bool=True, show_plot:bool=True,
                     fig=None, **kwargs_savefig):
        """
        The function saves the matplotlib figure.

        Parameters
        ----------
        fig_name : str
            name of the figure to use including figure extension.
        savefig : bool, optional
            Save the figure or not. The default is True.
        show_plot : bool, optional
            Show the plot or not. The default is True.
        fig : matplotlib figure instance, optional
            matplotlib figure instance. The default is None.
        **kwargs_savefig : matplotlib savefig kwargs
            Any other matplotlib savefig keywords arguments.

        Returns
        -------

        """
        return self._save_figure(fig_name, savefig=savefig, show_plot=show_plot,
                                 fig=fig, **kwargs_savefig)
