#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:05:16 2026

@author: badal.mondal
"""
from ._general_Hamiltonian_fns import _kpoints_from_path
import numpy as np

#%% ===========================================================================    
class _kp_H_30x30(_kpoints_from_path):
    # Default is return_eigen_val_vec_both=False => calculate only eigenvalues. 
    def __init__(self, kpath_list:str|list='L-G', nkpoints:int=41, 
                 return_eigen_val_vec_both:bool=False):
        #_kpoints_from_path.__init__(self)
        self.k_points = _kpoints_from_path._generate_k_points_from_path(kpath=kpath_list, nkpts=nkpoints)
        self.return_eigen_val_vec = return_eigen_val_vec_both
     
    @staticmethod
    def _generate_k_points_in_nm(a0, kpts_arr):
        # h_bar*h_bar*2*pi_*pi_/e_charge/e_mass*1e20 = 150.41206484194905
        # h^2.(2pi*k/(a*1e-10))^2/(2*m0) 
        k_sqr_const = 150.41206484194905/(a0*a0) 
        # 2*pi_/0.1 = 62.83185307179586
        # 2*pi/a*k_points [1/nm]
        k_const = 62.83185307179586/a0
        
        k_squared = k_sqr_const*np.sum(kpts_arr*kpts_arr, axis=1) # h^2.k^2/(2*m0) eV
        k_points =  k_const*kpts_arr # nm^-1
        return k_squared, k_points

    def _diagonalize_H_30x30(self, kp_params): 
        eigenvalvecs_k_path = {}
        for k_p_sec_name, kpts_arr in self.k_points.items():
            _k_squared, k_p_sec = self._generate_k_points_in_nm(kp_params['lattice_a0'], kpts_arr)
            eigenvalvecs = []
            for ii in range(len(k_p_sec)):
                 self.k = k_p_sec[ii]
                 self.k_sqr = _k_squared[ii]
                 eigenvalvecs.append(self._diagon_H_30x30(kp_params))
            # kx, ky, kz, e1, e2, e3, ...
            eigenvalvecs_k_path[k_p_sec_name] = np.concatenate((k_p_sec,eigenvalvecs), axis=1)
        return eigenvalvecs_k_path
    
    def _diagon_H_30x30(self, kp_params):
            return np.linalg.eigh(self._construct_H_30x30(kp_params)) if self.return_eigen_val_vec \
                    else np.linalg.eigvalsh(self._construct_H_30x30(kp_params))
        
    def _construct_H_30x30(self, kp_params):
        H_off_diag = np.zeros((30,30), dtype=np.complex128)
        H_diag = np.zeros((30,30), dtype=np.complex128)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Diagonal blocks
        H_diag[0:2, 0:2]     = self._H_diag_blocks(kp_params['EG_2u'], mat_dim=2) # H_2x2_2u
        H_diag[2:8, 2:8]     = self._H_diag_blocks(kp_params['EG_25u'], Edelta=kp_params['ED_25u'],
                                                   mat_dim=6) # H_6x6_25u
        H_diag[8:12, 8:12]   = self._H_diag_blocks(kp_params['EG_12'], mat_dim=4) # H_4x4_12
        H_diag[12:14, 12:14] = self._H_diag_blocks(kp_params['EG_1u'], mat_dim=2) # H_2x2_1u 
        H_diag[14:16, 14:16] = self._H_diag_blocks(kp_params['EG_1l'], mat_dim=2) # H_2x2_1l
        H_diag[16:22, 16:22] = self._H_diag_blocks(kp_params['EG_15'], Edelta=kp_params['ED_15'],
                                                   mat_dim=6) # H_6x6_15
        H_diag[22:24, 22:24] = self._H_diag_blocks(kp_params['EG_2l'], mat_dim=2) # H_2x2_2l
        H_diag[24:30, 24:30] = self._H_diag_blocks(kp_params['EG_25l'], Edelta=kp_params['ED_25l'],
                                                   mat_dim=6) # H_6x6_25l
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Off-diagonal block
        H_off_diag[0:2, 2:8] = self._Hk_2x6(mul_fact=kp_params['P_25u2u']) # H_2x6_p4
        H_off_diag[0:2, 16:22] = self._Hk_2x6(mul_fact=kp_params['P_152u']) # H_2x6_s2
        H_off_diag[0:2, 24:30] = self._Hk_2x6(mul_fact=kp_params['P_25l2u']) # H_2x6_p3
        H_off_diag[2:8, 8:12] = self._Hk_4x6(mul_fact=kp_params['P_25u12']).T # H_6x4_r2
        H_off_diag[2:8, 16:22] = self._Hk_6x6(mul_fact=kp_params['P_25u15']) # H_6x6_q2
        H_off_diag[2:8, 22:24] = self._Hk_2x6(mul_fact=kp_params['P_25u2l']).T # H_6x2_p2
        H_off_diag[8:12, 24:30] = self._Hk_4x6(mul_fact=kp_params['P_25l12']) # H_4x6_r1
        H_off_diag[12:14, 16:22] = self._Hk_2x6(mul_fact=kp_params['P_1u15']) # H_2x6_t1
        H_off_diag[14:16, 16:22] = self._Hk_2x6(mul_fact=kp_params['P_1l15']) # H_2x6_t2
        H_off_diag[16:22, 22:24] = self._Hk_2x6(mul_fact=kp_params['P_152l']).T # H_6x2_s1
        H_off_diag[16:22, 24:30] = self._Hk_6x6(mul_fact=kp_params['P_25l15']) + \
                                            self._H_G_SO(kp_params['ED_1525l']) # H_6x6_q1
        H_off_diag[22:24, 24:30] = self._Hk_2x6(mul_fact=kp_params['P_25l2l']) # H_2x6_p1
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Off-diagonal SO block
        H_off_diag[2:8, 24:30] = self._H_G_SO(kp_params['ED_25l25u'])
        
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