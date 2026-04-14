#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:05:16 2026

@author: badal.mondal
"""

import numpy as np

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

#%% ===========================================================================    
class _kp_H_30x30(_kpoints_from_path):
    # a0 in Angstrom
    # E_Gamma keys = (EG_2u, EG_25u, EG_12, EG_1u, EG_1l, EG_15, EG_2l, EG_25l) # eV
    # E_Delta keys = (ED_25u, ED_15, ED_25l, ED_25ul) # eV
    # Off_diag keys = (P1, P2, P3, P4, Q1, Q2, R1, R2, T1, T2) # eV.nm
    # Default is return_eigen_val_vec_both=False => calculate only eigenvalues. 
    def __init__(self, a0:float, E_Gamma:dict, E_Delta:dict, Off_diag:dict, 
                 kpath_list:str|list='L-G', nkpoints:int=41, 
                 return_eigen_val_vec_both:bool=False):
        #======================================================================
        # h_bar*h_bar*2*pi_*pi_/e_charge/e_mass*1e20 = 150.41206484194905
        # h^2.(2pi*k/(a*1e-10))^2/(2*m0) 
        k_sqr_const = 150.41206484194905/(a0*a0) 
        # 2*pi_/0.1 = 62.83185307179586
        # 2*pi/a*k_points [1/nm]
        k_const = 62.83185307179586/a0
        
        k_points = self._generate_k_points_from_path(kpath=kpath_list, nkpts=nkpoints) 
        self.k_squared, self.k_points = {}, {}
        for k_p_sec, kpts_arr in k_points.items():
            self.k_squared[k_p_sec] = k_sqr_const*np.sum(kpts_arr*kpts_arr, axis=1) # h^2.k^2/(2*m0) eV
            self.k_points[k_p_sec] =  k_const*kpts_arr
        #======================================================================
        self.E_Gamma = E_Gamma
        self.E_Delta = E_Delta
        self.Off_diag = Off_diag
        self.return_eigen_val_vec = return_eigen_val_vec_both

    def _diagonalize_H_30x30(self): 
        eigenvalvecs_k_path = {}
        for k_p_sec_name, k_p_sec in self.k_points.items():
            _k_squared = self.k_squared[k_p_sec_name]
            eigenvalvecs = []
            for ii in range(len(k_p_sec)):
                 self.k = k_p_sec[ii]
                 self.k_sqr = _k_squared[ii]
                 eigenvalvecs.append(self._diagon_H_30x30())
            # kx, ky, kz, e1, e2, e3, ...
            eigenvalvecs_k_path[k_p_sec_name] = np.concatenate((k_p_sec,eigenvalvecs), axis=1)
        return eigenvalvecs_k_path
    
    def _diagon_H_30x30(self):
            return np.linalg.eigh(self._construct_H_30x30()) if self.return_eigen_val_vec \
                    else np.linalg.eigvalsh(self._construct_H_30x30())
        
    def _construct_H_30x30(self):
        H_off_diag = np.zeros((30,30), dtype=np.complex128)
        H_diag = np.zeros((30,30), dtype=np.complex128)
        
        # Diagonal blocks
        H_diag[0:2, 0:2]     = self._H_diag_blocks(self.E_Gamma['EG_2u'], mat_dim=2) # H_2x2_2u
        H_diag[2:8, 2:8]     = self._H_diag_blocks(self.E_Gamma['EG_25u'], Edelta=self.E_Delta['ED_25u'], mat_dim=6) # H_6x6_25u
        H_diag[8:12, 8:12]   = self._H_diag_blocks(self.E_Gamma['EG_12'], mat_dim=4) # H_4x4_12
        H_diag[12:14, 12:14] = self._H_diag_blocks(self.E_Gamma['EG_1u'], mat_dim=2) # H_2x2_1u 
        H_diag[14:16, 14:16] = self._H_diag_blocks(self.E_Gamma['EG_1l'], mat_dim=2) # H_2x2_1l
        H_diag[16:22, 16:22] = self._H_diag_blocks(self.E_Gamma['EG_15'], Edelta=self.E_Delta['ED_15'], mat_dim=6) # H_6x6_15
        H_diag[22:24, 22:24] = self._H_diag_blocks(self.E_Gamma['EG_2l'], mat_dim=2) # H_2x2_2l
        H_diag[24:30, 24:30] = self._H_diag_blocks(self.E_Gamma['EG_25l'], Edelta=self.E_Delta['ED_25l'], mat_dim=6) # H_6x6_25l
        
        # Off-diagonal blocks
        H_off_diag[0:2, 2:8] = self._Hk_2x6(mul_fact=self.Off_diag['P4']) # H_2x6_p4
        H_off_diag[0:2, 24:30] = self._Hk_2x6(mul_fact=self.Off_diag['P3']) # H_2x6_p3
        H_off_diag[2:8, 8:12] = self._Hk_4x6(mul_fact=self.Off_diag['R2']).T # H_6x4_r2
        H_off_diag[2:8, 16:22] = self._Hk_6x6(mul_fact=self.Off_diag['Q2']) # H_6x6_q2
        H_off_diag[2:8, 22:24] = self._Hk_2x6(mul_fact=self.Off_diag['P2']).T # H_6x2_p2
        H_off_diag[2:8, 24:30] = self._H_G_SO(self.E_Delta['ED_25ul']) # H_6x6_SO_25ul
        H_off_diag[8:12, 24:30] = self._Hk_4x6(mul_fact=self.Off_diag['R1']) # H_4x6_r1
        H_off_diag[12:14, 16:22] = self._Hk_2x6(mul_fact=self.Off_diag['T1']) # H_2x6_t1
        H_off_diag[14:16, 16:22] = self._Hk_2x6(mul_fact=self.Off_diag['T2']) # H_2x6_t2
        H_off_diag[16:22, 24:30] = self._Hk_6x6(mul_fact=self.Off_diag['Q1']) # H_6x6_q1
        H_off_diag[22:24, 24:30] = self._Hk_2x6(mul_fact=self.Off_diag['P1']) # H_2x6_p1
        
        return H_off_diag + np.conjugate(H_off_diag).T + H_diag
    
    def _H_diag_blocks(self, Egamma, Edelta=0, mat_dim:int=2):
        # 2x2 diagonal Hamiltonian block
        diag_mat = np.zeros((mat_dim, mat_dim), dtype=np.complex128) 
        np.fill_diagonal(diag_mat, Egamma + self.k_sqr)
        return diag_mat + self._H_G_SO(Edelta) if mat_dim==6 else diag_mat
    
    def _Hk_6x6(self, mul_fact=1.0): 
        # 6X6 off-diagonal matrix
        XX = np.zeros((6,6), dtype=np.complex128)
        YY = np.array([[0, self.k[2], self.k[1]],
                       [self.k[2], 0, self.k[0]], 
                       [self.k[1], self.k[0], 0]])
        XX[0:3, 0:3] = YY
        XX[3:6, 3:6] = YY
        return XX*mul_fact
    
    def _Hk_4x6(self, mul_fact=1.0): 
        # 4x6 off-diagonal matrix
        sqrt3 = 1.7320508075688772
        XX = np.zeros((4,6), dtype=np.complex128)
        YY = np.array([[0, sqrt3*self.k[1], -sqrt3*self.k[2]],
                       [2.0*self.k[0], -self.k[1], -self.k[2]]])
        XX[0:2, 0:3] = YY
        XX[2:4, 3:6] = YY
        return XX*mul_fact
    
    def _Hk_2x6(self, mul_fact=1.0): 
        # 2x6 off-diagonal matrix
        XX = np.zeros((2,6), dtype=np.complex128)
        YY = np.array([self.k[0], self.k[1], self.k[2]])
        XX[0, 0:3] = YY
        XX[1, 3:6] = YY
        return XX*mul_fact
    
    @staticmethod
    def _H_G_SO(delta_G):
        # H_gamma_SO
        return delta_G/3.0 * np.array([[-1, -1j, 0,   0,  0,   1 ], 
                                       [1j, -1,  0,   0,  0,  -1j],
                                       [ 0,  0, -1,  -1,  1j,  0 ], 
                                       [ 0,  0, -1,  -1,  1j,  0 ],
                                       [ 0,  0, -1j, -1j, -1,  0 ],
                                       [ 1,  1j, 0,   0,   0, -1]])