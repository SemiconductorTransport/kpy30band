# Ref-1: Rideau et al., Phys. Rev. B 74, 195208 (2006)
# Ref-2: Song et. al. New J. Phys. 21, 073037 (2019) and 
# Ref-3: Song et. al. IEEE JOURNAL OF QUANTUM ELECTRONICS, 56, 7100208 (2020)

#==============================================================================
# NOTE: In nextnano all matrix elements are in a.u.. Here all matrix elements are
# in eV.nm unit. [conversion: X a.u = X*1.44/2 eV.nm]
# Mapping with nextnano:
#   E_1q=EG_2u,E_5d=EG_25u,E_3t= EG_12,E_1u=EG_1u
#   G_15,E_1c=EG_2l,E_1w=EG_1l
#   P_0=P_25l2l,P_1=P_25u2l,P_2=P_25l2u,P_3=P_25u2u,
#   P_4=P_1u15,P_5=P_1l15,P_prime_0=P_152l,P_prime_1=P_152u,
#   Q_0=P_25l15,Q_1=P_25u15,R_0=P_25l12,R_1=P_25u12 
#   delta_5v=ED_25l,delta_5c=ED_15,delta_5d=ED_25u,
#   delta_5v5c=ED_1525l,delta_5v5d=ED_25l25u 
#   l_5v=S_l_EG_25l,m_5v=S_m_EG_25l,n_5v=S_n_EG_25l 
#   l_5c=S_l_EG_15,m_5c=S_m_EG_15,n_5c=S_n_EG_15
#   l_5d=S_l_EG_25u,m_5d=S_m_EG_25u,n_5d=S_n_EG_25u
#   a_3t=S_a_12,b_3t=S_b_12,c_3t=S_c_12,d_3t=S_d_12
#   a_1c=S_a_EG_2l,a_1q=S_a_EG_2u,a_1w=S_a_EG_1l,a_1u=S_a_EG_1u 
#   l_5v5d=S_l_EG_25l25u,m_5v5d=S_m_EG_25l25u,n_5v5d=S_n_EG_25l25u
#   f_1u5d=S_f_EG_1u25u,f_1w5v=S_f_EG_1l25l,
#   f_1u5v=S_f_EG_1u25l, f_1w5d=S_f_EG_1l25u,
#   f_5c1c=S_f_EG_152l,f_5c1q=S_f_EG_152u
#   a_1c1q=S_a_EG_2l2u,a_1w1u=S_a_EG_1l1u,
#   g_3t1q=S_g_EG_122u,g_3t1c=S_a_EG_122l
#   h_5c3t=S_h_EG_1512
#==============================================================================

# Parameters nomenclature used here is based on References-2 and 3.
material_database = { 
    #==========================================================================
    # Binaries
    'test_sample': { # Taken from Ref-1
           # Lattice parameters
           'lattice_a0': 0, 'lattice_c0': 0,
           # E_Gamma (eV)
           'EG_1u': 0, 'EG_2u': 0, 'EG_25u': 0, 
           'EG_1l': 0, 'EG_2l': 0, 'EG_25l':0, 
           'EG_12': 0, 'EG_15': 0, 
           # E_Delta (eV)
           'ED_25u': 0, 'ED_15': 0, 'ED_25l': 0, 
           'ED_25l25u': 0, 'ED_1525l': 0,
           # Matrix element (eV.nm)
           'P_25l2l': 0, 'P_25u2l': 0, 'P_25l2u': 0, 'P_25u2u': 0, 
           'P_25l15': 0, 'P_25u15': 0, 'P_25l12': 0, 'P_25u12': 0, 
           'P_1u15': 0, 'P_1l15': 0, 'P_152l': 0, 'P_152u': 0,
           # Strain parameters
           'S_l_EG_25l': 0, 'S_m_EG_25l': 0, 'S_n_EG_25l': 0,
           'S_l_EG_15': 0, 'S_m_EG_15': 0, 'S_n_EG_15': 0,
           'S_l_EG_25u': 0, 'S_m_EG_25u': 0, 'S_n_EG_25u': 0,
           'S_a_12': 0, 'S_b_12': 5.47, 'S_c_12': 0, 'S_d_12': 0,
           'S_a_EG_2l': 0, 'S_a_EG_2u': 0, 'S_a_EG_1l': 0, 'S_a_EG_1u': 0,
           'S_l_EG_25l25u': 0, 'S_m_EG_25l25u': 0, 'S_n_EG_25l25u': 0,
           'S_f_EG_1u25u': 0, 'S_f_EG_1l25l': 0, 
           'S_f_EG_1u25l': 0, 'S_f_EG_1l25u': 0,
           'S_f_EG_152l': 0, 'S_f_EG_152u': 0,
           'S_a_EG_2l2u': 0,'S_a_EG_1l1u': 0,
           'S_g_EG_122u': 0, 'S_a_EG_122l': 0,
           'S_h_EG_1512': 0,
           # Elastic constants: C_ij parameters (GPa)
           'C_11': 0, 'C_12': 0, 'C_44': 0
	},
    'Si': { # Taken from Ref-1
           # Lattice parameters
           'lattice_a0': 5.387, 
           # E_Gamma (eV)
           'EG_1u': 8.4, 'EG_2u': 15.8, 'EG_25u': 11.7, 
           'EG_1l': -12.7, 'EG_2l': 4.15, 'EG_25l':0.00, 
           'EG_12': 8.54, 'EG_15': 3.335, 
           # E_Delta (eV)
           'ED_25u': 0.012, 'ED_15': 0.033, 'ED_25l': 0.044, 
           'ED_25l25u': 0.022, 'ED_1525l': 0.0,
           # Matrix element (eV.nm)
           'P_25l2l': 0.8784, 'P_25u2l': -0.00576, 'P_25l2u': 0.11232, 'P_25u2u': 1.026, 
           'P_25l15': 0.7689, 'P_25u15': -0.472, 'P_25l12': 0.3907, 'P_25u12': 0.6006, 
           'P_1u15': 0.8395, 'P_1l15': 0.2088, 'P_152l': 0.00, 'P_152u': 0.00,
           # Strain parameters
           'S_l_EG_25l': -2.7, 'S_m_EG_25l': 4.2, 'S_n_EG_25l': -7.379,
           'S_l_EG_15': 3.4, 'S_m_EG_15': -0.5, 'S_n_EG_15': -10.392,
           'S_l_EG_25u': -19.0, 'S_m_EG_25u': 8.0, 'S_n_EG_25u': -1.732,
           'S_a_12': 7.7, 'S_b_12': 5.47, 'S_c_12': 7.3, 'S_d_12': 3.65,
           'S_a_EG_2l': -9.0, 'S_a_EG_2u': 5.0, 'S_a_EG_1l': 10.0, 'S_a_EG_1u': 0.5,
           'S_l_EG_25l25u': -19.8, 'S_m_EG_25l25u': 3.9, 'S_n_EG_25l25u': 0.0,
           'S_f_EG_1u25u': 6.0, 'S_f_EG_1l25l': -5.0, 
           'S_f_EG_1u25l': -10.0, 'S_f_EG_1l25u': 0,
           'S_f_EG_152l': -19.0, 'S_f_EG_152u': -2.0,
           'S_a_EG_2l2u': 0.3,'S_a_EG_1l1u': -2,
           'S_g_EG_122u': -10.5,'S_a_EG_122l': -4.5,
           'S_h_EG_1512': 0,
           # Elastic constants: C_ij parameters (GPa)
           'C_11': 167.5, 'C_12': 65, 'C_44': 80.1
	},
    'Ge': { # Taken from Ref-2
           # Lattice parameters
           'lattice_a0': 5.658 , 
           # E_Gamma (eV)
           'EG_1u': 8.2064, 'EG_2u': 17.0426, 'EG_25u': 13.4041, 
           'EG_1l': -12.2519, 'EG_2l': 0.8140, 'EG_25l':0.00, 
           'EG_12': 8.5786, 'EG_15': 2.9900, 
           # E_Delta (eV)
           'ED_25u': 0.0793, 'ED_15': 0.2520, 'ED_25l': 0.2247, 
           'ED_25l25u': 0.22, 'ED_1525l': 0.0,
           # Matrix element (eV.nm)
           'P_25l2l': 0.8421, 'P_25u2l': 0.1781, 'P_25l2u': -0.0734, 'P_25u2u': 1.0543, 
           'P_25l15': 0.8114, 'P_25u15': -0.5334, 'P_25l12': 0.3757, 'P_25u12': 0.6820, 
           'P_1u15': 0.7994, 'P_1l15': -0.0384, 'P_152l': 0.00, 'P_152u': 0.00,
           # Strain parameters (Taken from Ref-3)
           'S_l_EG_25l': -3.8, 'S_m_EG_25l': 4.9, 'S_n_EG_25l': -9.527,
           'S_l_EG_15': 6.026, 'S_m_EG_15': 0.762, 'S_n_EG_15': -10.134,
           'S_l_EG_25u': -20.692, 'S_m_EG_25u': 9.119, 'S_n_EG_25u': 0.481,
           'S_a_12': 6.815, 'S_b_12': 6.798, 'S_c_12': 7.745, 'S_d_12': 4.858,
           'S_a_EG_2l': -7.181, 'S_a_EG_2u': 4.490, 'S_a_EG_1l': 14.171, 'S_a_EG_1u': -0.492,
           'S_l_EG_25l25u': -24.139, 'S_m_EG_25l25u': -0.124, 'S_n_EG_25l25u': -0.112,
           'S_f_EG_1u25u': 11.220, 'S_f_EG_1l25l': -7.666, 
           'S_f_EG_1u25l': 12.210, 'S_f_EG_1l25u': 0,
           'S_f_EG_152l': -22.242, 'S_f_EG_152u': 19.925,
           'S_a_EG_2l2u': -1.211,'S_a_EG_1l1u': -5.927,
           'S_g_EG_122u': -5.000,'S_a_EG_122l': -5.354,
           'S_h_EG_1512': 0,
           # Elastic constants: C_ij parameters (GPa)
           'C_11': 131.5, 'C_12': 49.4, 'C_44': 68.4
	},
    'Sn': { # Taken from Ref-2
           # Lattice parameters
           'lattice_a0': 6.489, 
           # E_Gamma (eV)
           'EG_1u': 5.473, 'EG_2u': 11.52, 'EG_25u': 8.546, 
           'EG_1l': -10.827, 'EG_2l': -0.376, 'EG_25l':0.00, 
           'EG_12': 7.593, 'EG_15': 2.194, 
           # E_Delta (eV)
           'ED_25u': 0.046, 'ED_15': 0.445, 'ED_25l': 0.652, 
           'ED_25l25u': 0.556, 'ED_1525l': 0.0,
           # Matrix element (eV.nm)
           'P_25l2l': 0.8425, 'P_25u2l': -0.0845, 'P_25l2u': 0.0014, 'P_25u2u': 0.8497, 
           'P_25l15':  0.6727, 'P_25u15': 0.4438, 'P_25l12': 0.4186, 'P_25u12': -0.4923, 
           'P_1u15': 0.6988, 'P_1l15': 0.3291, 'P_152l': 0.00, 'P_152u': 0.00,
           # Strain parameters (Taken from Ref-3)
           'S_l_EG_25l': 3.476, 'S_m_EG_25l': -0.047, 'S_n_EG_25l': 3.157,
           'S_l_EG_15': 40.493, 'S_m_EG_15': -39.341, 'S_n_EG_15': -1.437,
           'S_l_EG_25u': -6.156, 'S_m_EG_25u': 0.266, 'S_n_EG_25u': 0.535,
           'S_a_12': 0.457, 'S_b_12': -0.878, 'S_c_12': 0.192, 'S_d_12': -1.364,
           'S_a_EG_2l': -3.029, 'S_a_EG_2u': 16.445, 'S_a_EG_1l': -1.256, 'S_a_EG_1u': -16.983,
           'S_l_EG_25l25u': -42.944, 'S_m_EG_25l25u': 0.389, 'S_n_EG_25l25u': -2.369,
           'S_f_EG_1u25u': 6.17, 'S_f_EG_1l25l': -11.296, 
           'S_f_EG_1u25l': 7.293, 'S_f_EG_1l25u': 0,
           'S_f_EG_152l': 15.868, 'S_f_EG_152u': 44.263,
           'S_a_EG_2l2u': -2.592,'S_a_EG_1l1u': 8.543,
           'S_g_EG_122u': -1.721,'S_a_EG_122l': -9.579,
           'S_h_EG_1512': 0,
           # Elastic constants: C_ij parameters (GPa)
           'C_11': 69.0, 'C_12': 34.0, 'C_44': 30.0
	},
    #==========================================================================
    # Alloys
    'SiGe': { # Taken from Ref-1 (need to go through it again)
           # Lattice parameters
           'lattice_a0': 0, 
           # E_Gamma (eV)
           'EG_1u': 0.0000, 'EG_2u': 0.0000, 'EG_25u': 0.0000, 
           'EG_1l': 0.0000, 'EG_2l': 0.0000, 'EG_25l': 0.0000, 
           'EG_12': 0.0000, 'EG_15': 0.0000, 
           # E_Delta (eV)
           'ED_25u': 0.0000, 'ED_15': 0.0000, 'ED_25l': 0.0000, 
           'ED_25l25u': 0.0000, 'ED_1525l': 0.0000,
           # Matrix element (eV.nm)
           'P_25l2l': 0.0000, 'P_25u2l': 0.0000, 'P_25l2u': 0.0000, 'P_25u2u': 0.0000, 
           'P_25l15': 0.0000, 'P_25u15': 0.0000, 'P_25l12': 0.0000, 'P_25u12': 0.0000, 
           'P_1u15': 0.0000, 'P_1l15': 0.0000, 'P_152l': 0.0000, 'P_152u': 0.0000,
           # Strain parameters (Taken from Ref-3)
           'S_l_EG_25l': 0.0000, 'S_m_EG_25l': 0.0000, 'S_n_EG_25l': 0.0000,
           'S_l_EG_15': 0.0000, 'S_m_EG_15': 0.0000, 'S_n_EG_15': 0.0000,
           'S_l_EG_25u': 0.0000, 'S_m_EG_25u': 0.0000, 'S_n_EG_25u': 0.0000,
           'S_a_12': 0.0000, 'S_b_12': 0.0000, 'S_c_12': 0.0000, 'S_d_12': 0.0000,
           'S_a_EG_2l': 0.0000, 'S_a_EG_2u': 0.0000, 'S_a_EG_1l': 0.0000, 'S_a_EG_1u': 0.0000,
           'S_l_EG_25l25u': 0.0000, 'S_m_EG_25l25u': 0.0000, 'S_n_EG_25l25u': 0.0000,
           'S_f_EG_1u25u': 0.0000, 'S_f_EG_1l25l': 0.0000, 
           'S_f_EG_1u25l': 0.0000, 'S_f_EG_1l25u': 0.0000,
           'S_f_EG_152l': 0.0000, 'S_f_EG_152u': 0.0000,
           'S_a_EG_2l2u': 0.0000,'S_a_EG_1l1u': 0.0000,
           'S_g_EG_122u': 0.0000,'S_a_EG_122l': 0.0000,
           'S_h_EG_1512': 0.0000,
           'C_11': 0.0000, 'C_12': 0.0000, 'C_44': 0.0000
	},
    'GeSn': { # Taken from Ref-2
           # Lattice parameters
           'lattice_a0': 0, 
           # E_Gamma (eV)
           'EG_1u': 0.0000, 'EG_2u': 0.0000, 'EG_25u': 0.0000, 
           'EG_1l': 0.0000, 'EG_2l': 2.2767, 'EG_25l': 0.0000, 
           'EG_12': 0.0000, 'EG_15': 0.0000, 
           # E_Delta (eV)
           'ED_25u': 0.0000, 'ED_15': 0.0000, 'ED_25l': -4.9535, 
           'ED_25l25u': 0.0000, 'ED_1525l': 0.0000,
           # Matrix element (eV.nm)
           'P_25l2l': -0.3346, 'P_25u2l': 1.9328, 'P_25l2u': -1.2452, 'P_25u2u': -0.4220, 
           'P_25l15': -0.0055, 'P_25u15': 2.4647, 'P_25l12': 0.0089, 'P_25u12': 0.0000, 
           'P_1u15': -0.0441, 'P_1l15': 0.0000, 'P_152l': 0.0000, 'P_152u': 0.0000,
           # Strain parameters (Taken from Ref-3)
           'S_l_EG_25l': 0.0000, 'S_m_EG_25l': 0.0000, 'S_n_EG_25l': 0.0000,
           'S_l_EG_15': 0.0000, 'S_m_EG_15': 0.0000, 'S_n_EG_15': 0.0000,
           'S_l_EG_25u': 0.0000, 'S_m_EG_25u': 0.0000, 'S_n_EG_25u': 0.0000,
           'S_a_12': 0.0000, 'S_b_12': 0.0000, 'S_c_12': 0.0000, 'S_d_12': 0.0000,
           'S_a_EG_2l': 0.0000, 'S_a_EG_2u': 0.0000, 'S_a_EG_1l': 0.0000, 'S_a_EG_1u': 0.0000,
           'S_l_EG_25l25u': 0.0000, 'S_m_EG_25l25u': 0.0000, 'S_n_EG_25l25u': 0.0000,
           'S_f_EG_1u25u': 0.0000, 'S_f_EG_1l25l': 0.0000, 
           'S_f_EG_1u25l': 0.0000, 'S_f_EG_1l25u': 0.000,
           'S_f_EG_152l': 0.0000, 'S_f_EG_152u': 0.0000,
           'S_a_EG_2l2u': 0.0000,'S_a_EG_1l1u': 0.0000,
           'S_g_EG_122u': 0.0000,'S_a_EG_122l': 0.0000,
           'S_h_EG_1512': 0.0000,
           'C_11': 0.0000, 'C_12': 0.0000, 'C_44': 0.0000
	}
    }
