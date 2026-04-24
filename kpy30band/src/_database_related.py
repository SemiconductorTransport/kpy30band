#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:23:32 2026

@author: badal.mondal
"""

from .database import material_database

#%%# ==============================================================================
class _DataBase:
    def __init__(self):
        pass
    
    @classmethod
    def _print_database(cls, for_material=None):
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
        
        database_info = '''
Ref-1: Rideau et al., Phys. Rev. B 74, 195208 (2006)
Ref-2: Song et. al. New J. Phys. 21, 073037 (2019) and 
Ref-3: Song et. al. IEEE JOURNAL OF QUANTUM ELECTRONICS, 56, 7100208 (2020)
        '''
        database_header_info = """
The database contains the following information for each material. 
For ternary or higher-order materials, the parameter values represent the bowing coefficients.
        """
        if for_material is None:
            for_material = 'test_sample'
            print("""NOTE: Use 'print(DataBase().print_database(for_material=<material name> e.g. 'Si'))'
to print parameter values for the specific material""")
            print(f'{database_header_info}{database_info}\nParameters (parameter symbol : unit):')  

        params_db = material_database.get(for_material).copy()
        if params_db:
            database_header_info = f"""
NOTE: For ternary or higher-order materials, the parameter values represent the bowing coefficients.
The database contains the following information for {for_material}: 
            """
            print(database_header_info)
            for key, value in params_db.items():
                if key.startswith(('S_', 'ED_', 'EG_')):
                    unit = 'eV'
                elif key.startswith('P_'):
                    unit = 'eV.nm'
                elif key.startswith('C_'):
                    unit = 'GPa'
                elif key.startswith('lattice_'):    
                    unit = 'A'
                else:
                    unit= ''
                print(f'{key}: {value} {unit}')
        else:
            print(f'Error: "{for_material}" material does not exists in database yet. Contact developers.')
            
    @staticmethod
    def _print_materials_available_in_database():
        database_footer_info = """
NOTE: For ternary or higher-order materials, the parameter values represent the bowing coefficients.
        """
        avail_material_list = ', '.join(list(material_database.keys()))
        print(f'* Available materials name: {avail_material_list}{database_footer_info}')
        
        
        