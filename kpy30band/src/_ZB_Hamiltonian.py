#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:05:16 2026

@author: badal.mondal
"""
from ._general_Hamiltonian_fns import _kpoints_from_point_path
import numpy as np

#%% ===========================================================================    
class _kp_H_30x30(_kpoints_from_point_path):
    # Default is return_eigen_val_vec_both=False => calculate only eigenvalues. 
    def __init__(self, kpoints_list:str|list='L-G', nkpoints:int=41,
                 return_eigen_val_vec_both:bool=False):
        #_kpoints_from_path.__init__(self)
        self.k_points = _kpoints_from_point_path._generate_k_points_from_point_path(kpath=kpoints_list, 
                                                                                    nkpts=nkpoints)
        self.return_eigen_val_vec = return_eigen_val_vec_both
     
    @staticmethod
    def _generate_k_points_in_nm(a0, kpts_arr):
        # h_bar*h_bar*2*pi_*pi_/e_charge/e_mass*1e20 = 150.41206484194905
        # h^2.(2pi*k/(a*1e-10))^2/(2*m0) 
        k_sqr_const = 150.41206484194905/(a0*a0) 
        # 2*pi_/0.1 = 62.83185307179586
        # 2*pi/a*k_points [1/nm]
        k_const = 62.83185307179586/a0
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
        if self.return_eigen_val_vec:
            return np.linalg.eigh(self._construct_H_30x30(kp_params))  
        return np.linalg.eigvalsh(self._construct_H_30x30(kp_params))
                    
    def _construct_H_30x30(self, kp_params):
        if self.apply_strain_:
            return self._construct_unstrained_H_30x30(kp_params) + \
                            self._construct_strained_H_30x30(kp_params)
        return self._construct_unstrained_H_30x30(kp_params)
    
    #**************************************************************************
    def _construct_unstrained_H_30x30(self, kp_params):
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
        H_off_diag[0:2, 2:8]     = self._Hk_2x6(mul_fact=kp_params['P_25u2u']) # H_2x6_p4
        H_off_diag[0:2, 16:22]   = self._Hk_2x6(mul_fact=kp_params['P_152u']) # H_2x6_s2
        H_off_diag[0:2, 24:30]   = self._Hk_2x6(mul_fact=kp_params['P_25l2u']) # H_2x6_p3
        H_off_diag[2:8, 8:12]    = self._Hk_4x6(mul_fact=kp_params['P_25u12']).T # H_6x4_r2
        H_off_diag[2:8, 16:22]   = self._Hk_6x6(mul_fact=kp_params['P_25u15']) # H_6x6_q2
        H_off_diag[2:8, 22:24]   = self._Hk_2x6(mul_fact=kp_params['P_25u2l']).T # H_6x2_p2
        H_off_diag[8:12, 24:30]  = self._Hk_4x6(mul_fact=kp_params['P_25l12']) # H_4x6_r1
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
    #**************************************************************************
    def _construct_strained_H_30x30(self, kp_params):
        W_off_diag = np.zeros((30,30), dtype=np.complex128)
        W_diag = np.zeros((30,30), dtype=np.complex128)
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Diagonal blocks
        W_diag[0:2, 0:2]     = self._W_diag_blocks_2x2(kp_params['S_a_EG_2u']) 
        W_diag[2:8, 2:8]     = self._W_diag_blocks_6x6(kp_params['S_l_EG_25u'],
                                                       kp_params['S_m_EG_25u'],
                                                       kp_params['S_n_EG_25u'])
        W_diag[8:12, 8:12]   = self._W_diag_blocks_4x4(kp_params['S_a_12'],
                                                       kp_params['S_b_12'],
                                                       kp_params['S_c_12'],
                                                       kp_params['S_d_12'])
        W_diag[12:14, 12:14] = self._W_diag_blocks_2x2(kp_params['S_a_EG_1u']) 
        W_diag[14:16, 14:16] = self._W_diag_blocks_2x2(kp_params['S_a_EG_1l'])
        W_diag[16:22, 16:22] = self._W_diag_blocks_6x6(kp_params['S_l_EG_15'],
                                                       kp_params['S_m_EG_15'],
                                                       kp_params['S_n_EG_15'])
        W_diag[22:24, 22:24] = self._W_diag_blocks_2x2(kp_params['S_a_EG_2l']) 
        W_diag[24:30, 24:30] = self._W_diag_blocks_6x6(kp_params['S_l_EG_25l'],
                                                       kp_params['S_m_EG_25l'],
                                                       kp_params['S_n_EG_25l'])
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Off-diagonal k-independent blocks
        W_off_diag[0:2, 8:12]    = self._W_off_diag_blocks_4x2(kp_params['S_g_EG_122u']).T
        W_off_diag[0:2, 16:22]   = self._W_off_diag_blocks_2x6(kp_params['S_f_EG_152u'])
        W_off_diag[0:2, 22:24]   = self._W_diag_blocks_2x2(kp_params['S_a_EG_2l2u'])
        W_off_diag[2:8, 14:16]   = self._W_off_diag_blocks_2x6(kp_params['S_f_EG_1u25u']).T
        W_off_diag[2:8, 16:18]   = self._W_off_diag_blocks_2x6(kp_params['S_f_EG_1l25u']).T
        W_off_diag[2:8, 24:30]   = self._W_diag_blocks_6x6(kp_params['S_l_EG_25l25u'],
                                                           kp_params['S_m_EG_25l25u'],
                                                           kp_params['S_n_EG_25l25u'])
        W_off_diag[8:12, 16:22]  = self._W_off_diag_blocks_4x6(kp_params['S_h_EG_1512']) 
        W_off_diag[8:12, 22:24]  = self._W_off_diag_blocks_4x2(kp_params['S_a_EG_122l'])
        W_off_diag[12:14, 14:16] = self._W_diag_blocks_2x2(kp_params['S_a_EG_1l1u'])
        W_off_diag[12:14, 24:30] = self._W_off_diag_blocks_2x6(kp_params['S_f_EG_1u25l'])
        W_off_diag[14:16, 24:30] = self._W_off_diag_blocks_2x6(kp_params['S_f_EG_1l25l'])
        W_off_diag[16:22, 22:24] = self._W_off_diag_blocks_2x6(kp_params['S_f_EG_152l']).T
        
        # Off-diagonal block k-dependent blocks
        W_6x6, W_4x6, W_2x6 = self._Wk_off_diagonal()
        
        W_off_diag[0:2, 2:8]     = kp_params['P_25u2u'] * W_2x6
        W_off_diag[0:2, 24:30]   = kp_params['P_25l2u'] * W_2x6
        W_off_diag[2:8, 8:12]    = kp_params['P_25u12'] * W_4x6.T 
        W_off_diag[2:8, 16:22]   = kp_params['P_25u15'] * W_6x6
        W_off_diag[2:8, 22:24]   = kp_params['P_25u2l'] * W_2x6.T 
        W_off_diag[8:12, 24:30]  = kp_params['P_25l12'] * W_4x6
        W_off_diag[12:14, 16:22] = kp_params['P_1u15']  * W_2x6
        W_off_diag[14:16, 16:22] = kp_params['P_1l15']  * W_2x6
        W_off_diag[16:22, 24:30] = kp_params['P_25l15'] * W_6x6 
        W_off_diag[22:24, 24:30] = kp_params['P_25l2l'] * W_2x6 

        return W_off_diag + np.conjugate(W_off_diag).T + W_diag
    
    def _W_diag_blocks_2x2(self, a_Gamma):
        # 2x2 diagonal Hamiltonian block
        diag_mat = np.zeros((2, 2), dtype=np.complex128) 
        np.fill_diagonal(diag_mat, self.strain_tensor[0]+
                                   self.strain_tensor[1]+
                                   self.strain_tensor[2])
        return a_Gamma * diag_mat
    
    def _W_diag_blocks_6x6(self, l, m, n):
        # 6x6 diagonal Hamiltonian block
        # strain_tensor = [e_xx, e_yy, e_zz, e_yz, e_xz, e_xy]
        diag_mat = np.zeros((6, 6), dtype=np.complex128) 
        W_diag_blocks_3x3 = np.array([[l*self.strain_tensor[0]+
                                      m*(self.strain_tensor[1]+self.strain_tensor[2]), 
                                      n*self.strain_tensor[5], n*self.strain_tensor[4]],
                                     [n*self.strain_tensor[5], 
                                      l*self.strain_tensor[1]+
                                      m*(self.strain_tensor[0]+self.strain_tensor[2]),
                                      n*self.strain_tensor[3]],
                                     [n*self.strain_tensor[4], n*self.strain_tensor[3],
                                      l*self.strain_tensor[2]+
                                      m*(self.strain_tensor[0]+self.strain_tensor[1])]])
        diag_mat[0:3, 0:3] = W_diag_blocks_3x3
        diag_mat[3:6, 3:6] = W_diag_blocks_3x3
        return diag_mat
    
    def _W_diag_blocks_4x4(self, a, b, c, d):
        # 4x4 diagonal Hamiltonian block
        # strain_tensor = [e_xx, e_yy, e_zz, e_yz, e_xz, e_xy]
        A = 6.0*(b-d); B = 3.0*(a+b-2.0*c); C = 2.0*(2.0*a-4.0*c+b+d)
        D = 5.0*b-2.0*c-4.0*d+a; E = 1.7320508075688772*(2.0*c-2.0*d-a+b)
        
        diag_mat = np.zeros((4, 4), dtype=np.complex128) 
        e_11 = A*self.strain_tensor[0] + B*(self.strain_tensor[1]+self.strain_tensor[2])
        e_22 = C*self.strain_tensor[0] + D*(self.strain_tensor[1]+self.strain_tensor[2])
        e_12 = E*(self.strain_tensor[1]-self.strain_tensor[2])
        W_diag_blocks_2x2 = np.array([[e_11, e_12], [e_12, e_22]])
        diag_mat[0:2, 0:2] = W_diag_blocks_2x2
        diag_mat[2:4, 2:4] = W_diag_blocks_2x2
        return diag_mat
    
    def _W_off_diag_blocks_4x2(self, g_Gamma):
        # 4x2 off-diagonal Hamiltonian block
        # strain_tensor = [e_xx, e_yy, e_zz, e_yz, e_xz, e_xy]
        e_11 = 1.7320508075688772*(self.strain_tensor[1]-self.strain_tensor[2])
        e_12 = 2.0*self.strain_tensor[0]-self.strain_tensor[1]-self.strain_tensor[2]
        return g_Gamma * np.array([[e_11, 0], [e_12, 0], [0, e_11], [0, e_12]])
    
    def _W_off_diag_blocks_2x6(self, f_Gamma):
        # 2x6 off-diagonal Hamiltonian block
        # strain_tensor = [e_xx, e_yy, e_zz, e_yz, e_xz, e_xy]
        return f_Gamma * np.array([[self.strain_tensor[3], self.strain_tensor[4], 
                                    self.strain_tensor[5], 0, 0, 0], 
                                   [0, 0, 0, self.strain_tensor[3], 
                                    self.strain_tensor[4], self.strain_tensor[5]]])
    
    def _W_off_diag_blocks_4x6(self, h_Gamma):
        # 4x6 off diagonal Hamiltonian block
        
        diag_mat = np.zeros((4, 6), dtype=np.complex128) 
        sqrt3 = 1.7320508075688772
        W_2x3 = np.array([[0, sqrt3*self.strain_tensor[4], -1*sqrt3*self.strain_tensor[5]], 
                          [2.0*self.strain_tensor[3], -1*self.strain_tensor[4], -1*self.strain_tensor[5]]])
        diag_mat[0:2, 0:3] = W_2x3
        diag_mat[2:4, 3:6] = W_2x3
        return h_Gamma * diag_mat
    
    def _Wk_off_diagonal(self): 
        sqrt3 = 1.7320508075688772
        W_6x6 = np.zeros((6,6), dtype=np.complex128)
        W_4x6 = np.zeros((4,6), dtype=np.complex128)
        W_2x6 = np.zeros((2,6), dtype=np.complex128)
        
        W_3x3 = np.zeros((3,3), dtype=np.complex128)
        W_2x3 = np.zeros((2,3), dtype=np.complex128)
        W_1x3 = np.zeros((1,3), dtype=np.complex128)
        # strain_tensor = [e_xx, e_yy, e_zz, e_yz, e_xz, e_xy]
        # mapping w.r.t position in the strain_tensor
        e_iz_id = [4, 3, 2] # exz, eyz, ezz
        e_iy_id = [5, 1, 3] # e_xy, eyy, ezy
        e_ix_id = [0, 5, 4] # e_xx, eyx, ezx
        
        for i in range(3): # loop over x,y,z
            W_3x3 += np.array([[0, self.strain_tensor[e_iz_id[i]]*self.k[i], 
                                self.strain_tensor[e_iy_id[i]]*self.k[i]],
                               [self.strain_tensor[e_iz_id[i]]*self.k[i], 0, 
                                self.strain_tensor[e_ix_id[i]]*self.k[i]], 
                               [self.strain_tensor[e_iy_id[i]]*self.k[i],  
                                self.strain_tensor[e_ix_id[i]]*self.k[i], 0]])
            W_2x3 += np.array([[0, sqrt3*self.strain_tensor[e_iy_id[i]]*self.k[i], 
                                -sqrt3*self.strain_tensor[e_iz_id[i]]*self.k[i]],
                               [2.0*self.strain_tensor[e_ix_id[i]]*self.k[i],  
                                -1*self.strain_tensor[e_iy_id[i]]*self.k[i],
                                -1*self.strain_tensor[e_iz_id[i]]*self.k[i]]])
            W_1x3 += np.array([[self.strain_tensor[e_ix_id[i]]*self.k[i], 
                                self.strain_tensor[e_iy_id[i]]*self.k[i],
                                self.strain_tensor[e_iz_id[i]]*self.k[i]]])
            
        W_6x6[0:3, 0:3] = W_3x3
        W_6x6[3:6, 3:6] = W_3x3
        
        W_4x6[0:2, 0:3] = W_2x3
        W_4x6[2:4, 3:6] = W_2x3
        
        W_2x6[0:1, 0:3] = W_1x3
        W_2x6[1:2, 3:6] = W_1x3
        
        return -1.0*W_6x6, -1.0*W_4x6, -1.0*W_2x6