# Ref-1: Rideau et al., Phys. Rev. B 74, 195208 (2006)
# Ref-2: Song et. al. New J. Phys. 21, 073037 (2019) 
# Ref-3: Song et. al. IEEE JOURNAL OF QUANTUM ELECTRONICS, 56, 7100208 (2020)

#==============================================================================
# NOTE: In nextnano all matrix elements are in a.u.. Here all matrix elements are
# in eV.nm unit. [conversion: X a.u = X*hbar^2/(m0*a0) = X*1.44/2 eV.nm]
# Mapping with nextnano:
#   E_1q=EG_1q,E_5d=EG_5d,E_3t= EG_3,E_1u=EG_1u
#   EG_5c=EG_5c,E_1c=EG_1,E_1w=EG_1v,E_5v=EG_5
#   P_0=P_P,P_1=P_Pd,P_2=P_P2,P_3=P_P2d,
#   P_4=P_Pu,P_5=P_Ps,P_prime_0=P_Pprime,P_prime_1=P_P2prime,
#   Q_0=P_Px,Q_1=P_Pxd,R_0=P_P3,R_1=P_P3d 
#   delta_5v=ED_so,delta_5c=ED_c,delta_5d=ED_d,
#   delta_5v5c=ED_prime,delta_5v5d=ED_dso 
#   l_5v=S_l_EG_5,m_5v=S_m_EG_5,n_5v=S_n_EG_5 
#   l_5c=S_l_EG_5c,m_5c=S_m_EG_5c,n_5c=S_n_EG_5c
#   l_5d=S_l_EG_5d,m_5d=S_m_EG_5d,n_5d=S_n_EG_5d
#   a_3t=S_a_3,b_3t=S_b_3,c_3t=S_c_3,d_3t=S_d_3
#   a_1c=S_a_EG_1,a_1q=S_a_EG_1q,a_1w=S_a_EG_1v,a_1u=S_a_EG_1u 
#   l_5v5d=S_l_EG_55d,m_5v5d=S_m_EG_55d,n_5v5d=S_n_EG_55d
#   f_1u5d=S_f_EG_1u5d,f_1w5v=S_f_EG_1v5,
#   f_1u5v=S_f_EG_1u5, f_1w5d=S_f_EG_1v5d,
#   f_5c1c=S_f_EG_5c1,f_5c1q=S_f_EG_5c1q
#   a_1c1q=S_a_EG_11q,a_1w1u=S_a_EG_1v1u,
#   g_3t1q=S_g_EG_31q,g_3t1c=S_g_EG_31
#   h_5c3t=S_h_EG_5c3
#==============================================================================

# Parameters nomenclature used here is based on References-2 and 3.
material_database = { 
    #==========================================================================
    # Binaries
    'test_sample': { # 
                    'comment': 'Test sample',
                    # Lattice parameters
                    'lattice_a0': 0, 'lattice_c0': 0,
                    # E_Gamma (eV)
                    'EG_1u': 0, 'EG_1q': 0, 'EG_5d': 0, 
                    'EG_1v': 0, 'EG_1': 0, 'EG_5':0, 
                    'EG_3': 0, 'EG_5c': 0, 
                    # E_Delta (eV)
                    'ED_d': 0, 'ED_c': 0, 'ED_so': 0, 
                    'ED_dso': 0, 'ED_prime': 0, 'ED_cd_prime': 0,
                    # Matrix element (eV.nm)
                    'P_P': 0, 'P_Pd': 0, 'P_P2': 0, 'P_P2d': 0, 
                    'P_Px': 0, 'P_Pxd': 0, 'P_P3': 0, 'P_P3d': 0, 
                    'P_Pu': 0, 'P_Ps': 0, 'P_P_prime': 0, 'P_P2_prime': 0,
                    # Strain parameters
                    'S_l_EG_5': 0, 'S_m_EG_5': 0, 'S_n_EG_5': 0,
                    'S_l_EG_5c': 0, 'S_m_EG_5c': 0, 'S_n_EG_5c': 0,
                    'S_l_EG_5d': 0, 'S_m_EG_5d': 0, 'S_n_EG_5d': 0,
                    'S_a_3': 0, 'S_b_3': 5.47, 'S_c_3': 0, 'S_d_3': 0,
                    'S_a_EG_1': 0, 'S_a_EG_1q': 0, 'S_a_EG_1v': 0, 'S_a_EG_1u': 0,
                    'S_l_EG_55d': 0, 'S_m_EG_55d': 0, 'S_n_EG_55d': 0,
                    'S_f_EG_1u5d': 0, 'S_f_EG_1v5': 0, 
                    'S_f_EG_1u5': 0, 'S_f_EG_1v5d': 0,
                    'S_f_EG_5c1': 0, 'S_f_EG_5c1q': 0,
                    'S_a_EG_11q': 0,'S_a_EG_1v1u': 0,
                    'S_g_EG_31q': 0, 'S_g_EG_31': 0,
                    'S_h_EG_5c3': 0,
                    # Elastic constants: C_ij parameters (GPa)
                    'C_11': 0, 'C_12': 0, 'C_44': 0
                    }, 
    'Si': { # Taken from Ref-1
           'comment': 'Taken from Rideau et al., Phys. Rev. B 74, 195208 (2006)',
           # Lattice parameters
           'lattice_a0': 5.387, 
           # E_Gamma (eV)
           'EG_1u': 8.4, 'EG_1q': 15.8, 'EG_5d': 11.7, 
           'EG_1v': -12.7, 'EG_1': 4.15, 'EG_5':0.00, 
           'EG_3': 8.54, 'EG_5c': 3.335, 
           # E_Delta (eV)
           'ED_d': 0.012, 'ED_c': 0.033, 'ED_so': 0.044, 
           'ED_dso': 0.022, 'ED_prime': 0.0, 'ED_cd_prime': 0.0000,
           # Matrix element (eV.nm)
           'P_P': 0.8784, 'P_Pd': -0.0058, 'P_P2': 0.1123, 'P_P2d': 1.0260, 
           'P_Px': 0.7689, 'P_Pxd': -0.4720, 'P_P3': 0.3907, 'P_P3d': 0.6006, 
           'P_Pu': 0.8395, 'P_Ps': 0.2088, 'P_P_prime': 0.00, 'P_P2_prime': 0.00,
           # Strain parameters
           'S_l_EG_5': -2.7, 'S_m_EG_5': 4.2, 'S_n_EG_5': -7.379,
           'S_l_EG_5c': 3.4, 'S_m_EG_5c': -0.5, 'S_n_EG_5c': -10.392,
           'S_l_EG_5d': -19.0, 'S_m_EG_5d': 8.0, 'S_n_EG_5d': -1.732,
           'S_a_3': 7.7, 'S_b_3': 5.47, 'S_c_3': 7.3, 'S_d_3': 3.65,
           'S_a_EG_1': -9.0, 'S_a_EG_1q': 5.0, 'S_a_EG_1v': 10.0, 'S_a_EG_1u': 0.5,
           'S_l_EG_55d': -19.8, 'S_m_EG_55d': 3.9, 'S_n_EG_55d': 0.0,
           'S_f_EG_1u5d': 6.0, 'S_f_EG_1v5': -5.0, 
           'S_f_EG_1u5': -10.0, 'S_f_EG_1v5d': 0,
           'S_f_EG_5c1': -19.0, 'S_f_EG_5c1q': -2.0,
           'S_a_EG_11q': 0.3,'S_a_EG_1v1u': -2,
           'S_g_EG_31q': -10.5,'S_g_EG_31': -4.5,
           'S_h_EG_5c3': 0,
           # Elastic constants: C_ij parameters (GPa)
           'C_11': 168.3, 'C_12': 66.8, 'C_44': 79.9	
           },
    # 'Ge2' is same as 'Ge'. Because we have two sets of parameters for Ge.
    'Ge2': { 
           'comment': 'Taken from Rideau et al., Phys. Rev. B 74, 195208 (2006)',
           # Lattice parameters
           'lattice_a0': 5.583 , 
           # E_Gamma (eV)
           'EG_1u': 6.8, 'EG_1q': 14.0, 'EG_5d': 11.36, 
           'EG_1v': -12.88, 'EG_1': 0.89, 'EG_5': 0.00, 
           'EG_3': 10.3, 'EG_5c': 3.113, 
           # E_Delta (eV)
           'ED_d': 0.042, 'ED_c': 0.19, 'ED_so': 0.296, 
           'ED_dso': 0.22, 'ED_prime': 0.0, 'ED_cd_prime': 0.0000,
           # Matrix element (eV.nm)
           'P_P': 0.8539, 'P_Pd': 0.0144, 'P_P2': 0.1065, 'P_P2d': 1.0071, 
           'P_Px': 0.7738, 'P_Pxd': -0.5477, 'P_P3': 0.4544, 'P_P3d': 0.5915, 
           'P_Pu': 0.7929, 'P_Ps': 0.2664 , 'P_P_prime': 0.00, 'P_P2_prime': 0.00,
           # Strain parameters (Taken from Ref-3)
           'S_l_EG_5': -3.8, 'S_m_EG_5': 4.9, 'S_n_EG_5': -9.527,
           'S_l_EG_5c': 6.026, 'S_m_EG_5c': 0.762, 'S_n_EG_5c': -10.134,
           'S_l_EG_5d': -20.692, 'S_m_EG_5d': 9.119, 'S_n_EG_5d': 0.481,
           'S_a_3': 6.815, 'S_b_3': 6.798, 'S_c_3': 7.745, 'S_d_3': 4.858,
           'S_a_EG_1': -7.181, 'S_a_EG_1q': 4.490, 'S_a_EG_1v': 14.171, 'S_a_EG_1u': -0.492,
           'S_l_EG_55d': -24.139, 'S_m_EG_55d': -0.124, 'S_n_EG_55d': -0.112,
           'S_f_EG_1u5d': 11.220, 'S_f_EG_1v5': -7.666, 
           'S_f_EG_1u5': -12.210, 'S_f_EG_1v5d': 0,
           'S_f_EG_5c1': -22.242, 'S_f_EG_5c1q': 19.925,
           'S_a_EG_11q': -1.211,'S_a_EG_1v1u': -5.927,
           'S_g_EG_31q': -5.000,'S_g_EG_31': -5.354,
           'S_h_EG_5c3': 0,
           # Elastic constants: C_ij parameters (GPa)
           'C_11': 132.8, 'C_12': 46.8, 'C_44': 66.57
	},
    'Ge': { 
           'comment': 'Taken from Song et. al. New J. Phys. 21, 073037 (2019)',
           # Lattice parameters
           'lattice_a0': 5.658 , 
           # E_Gamma (eV)
           'EG_1u': 8.2064, 'EG_1q': 17.0426, 'EG_5d': 13.4041, 
           'EG_1v': -12.2519, 'EG_1': 0.8140, 'EG_5':0.00, 
           'EG_3': 8.5786, 'EG_5c': 2.9900, 
           # E_Delta (eV)
           'ED_d': 0.0793, 'ED_c': 0.2520, 'ED_so': 0.2247, 
           'ED_dso': 0.22, 'ED_prime': 0.0, 'ED_cd_prime': 0.0000,
           # Matrix element (eV.nm)
           'P_P': 0.8421, 'P_Pd': 0.1781, 'P_P2': -0.0734, 'P_P2d': 1.0543, 
           'P_Px': 0.8114, 'P_Pxd': -0.5334, 'P_P3': 0.3757, 'P_P3d': 0.6820, 
           'P_Pu': 0.7994, 'P_Ps': -0.0384, 'P_P_prime': 0.00, 'P_P2_prime': 0.00,
           # Strain parameters (Taken from Ref-3)
           'S_l_EG_5': -3.8, 'S_m_EG_5': 4.9, 'S_n_EG_5': -9.527,
           'S_l_EG_5c': 6.026, 'S_m_EG_5c': 0.762, 'S_n_EG_5c': -10.134,
           'S_l_EG_5d': -20.692, 'S_m_EG_5d': 9.119, 'S_n_EG_5d': 0.481,
           'S_a_3': 6.815, 'S_b_3': 6.798, 'S_c_3': 7.745, 'S_d_3': 4.858,
           'S_a_EG_1': -7.181, 'S_a_EG_1q': 4.490, 'S_a_EG_1v': 14.171, 'S_a_EG_1u': -0.492,
           'S_l_EG_55d': -24.139, 'S_m_EG_55d': -0.124, 'S_n_EG_55d': -0.112,
           'S_f_EG_1u5d': 11.220, 'S_f_EG_1v5': -7.666, 
           'S_f_EG_1u5': -12.210, 'S_f_EG_1v5d': 0,
           'S_f_EG_5c1': -22.242, 'S_f_EG_5c1q': 19.925,
           'S_a_EG_11q': -1.211,'S_a_EG_1v1u': -5.927,
           'S_g_EG_31q': -5.000,'S_g_EG_31': -5.354,
           'S_h_EG_5c3': 0,
           # Elastic constants: C_ij parameters (GPa)
           'C_11': 124.0, 'C_12': 41.3, 'C_44': 68.3
	},
    'Sn': { 
           'comment': 'Taken from Song et. al. New J. Phys. 21, 073037 (2019)',
           # Lattice parameters
           'lattice_a0': 6.4892, 
           # E_Gamma (eV)
           'EG_1u': 5.473, 'EG_1q': 11.52, 'EG_5d': 8.546, 
           'EG_1v': -10.827, 'EG_1': -0.376, 'EG_5':0.00, 
           'EG_3': 7.593, 'EG_5c': 2.194, 
           # E_Delta (eV)
           'ED_d': 0.046, 'ED_c': 0.445, 'ED_so': 0.652, 
           'ED_dso': 0.556, 'ED_prime': 0.0, 'ED_cd_prime': 0.0000,
           # Matrix element (eV.nm)
           'P_P': 0.8425, 'P_Pd': -0.0845, 'P_P2': 0.0014, 'P_P2d': 0.8497, 
           'P_Px':  0.6727, 'P_Pxd': 0.4438, 'P_P3': 0.4186, 'P_P3d': -0.4923, 
           'P_Pu': 0.6988, 'P_Ps': 0.3291, 'P_P_prime': 0.00, 'P_P2_prime': 0.00,
           # Strain parameters (Taken from Ref-3)
           'S_l_EG_5': 3.476, 'S_m_EG_5': -0.047, 'S_n_EG_5': 3.157,
           'S_l_EG_5c': 40.493, 'S_m_EG_5c': -39.341, 'S_n_EG_5c': -1.437,
           'S_l_EG_5d': -6.156, 'S_m_EG_5d': 0.266, 'S_n_EG_5d': 0.535,
           'S_a_3': 0.457, 'S_b_3': -0.878, 'S_c_3': 0.192, 'S_d_3': -1.364,
           'S_a_EG_1': -3.029, 'S_a_EG_1q': 16.445, 'S_a_EG_1v': -1.256, 'S_a_EG_1u': -16.983,
           'S_l_EG_55d': -42.944, 'S_m_EG_55d': 0.389, 'S_n_EG_55d': -2.369,
           'S_f_EG_1u5d': 6.17, 'S_f_EG_1v5': -11.296, 
           'S_f_EG_1u5': 7.293, 'S_f_EG_1v5d': 0,
           'S_f_EG_5c1': 15.868, 'S_f_EG_5c1q': 44.263,
           'S_a_EG_11q': -2.592,'S_a_EG_1v1u': 8.543,
           'S_g_EG_31q': -1.721,'S_g_EG_31': -9.579,
           'S_h_EG_5c3': 0,
           # Elastic constants: C_ij parameters (GPa)
           'C_11': 69.0, 'C_12': 29.3, 'C_44': 36.2
	},
    #==========================================================================
    # Alloys
    'SiGe2': { 
           'comment': 'Taken from Rideau et al., Phys. Rev. B 74, 195208 (2006)',
           # Lattice parameters
           'lattice_a0': 0.0532, 
           # E_Gamma (eV)
           'EG_1u': 0.0000, 'EG_1q': 0.0000, 'EG_5d': 0.0000, 
           'EG_1v': 0.0000, 'EG_1': 0.0000, 'EG_5': 0.0000, 
           'EG_3': 0.0000, 'EG_5c': 0.0000, 
           # E_Delta (eV)
           'ED_d': 0.0000, 'ED_c': 0.0000, 'ED_so': 0.052, 
           'ED_dso': 0.0000, 'ED_prime': -0.04, 'ED_cd_prime': 0.0000,
           # Matrix element (eV.nm)
           'P_P': 0.0000, 'P_Pd': -0.0360, 'P_P2': 0.0000, 'P_P2d': 0.0000, 
           'P_Px': 0.0000, 'P_Pxd': 0.0000, 'P_P3': 0.0000, 'P_P3d': 0.0000, 
           'P_Pu': -0.0288, 'P_Ps': 0.0000, 'P_P_prime': -0.1j, 'P_P2_prime': 0.3j,
           # Strain parameters (Taken from Ref-3)
           'S_l_EG_5': 0.0000, 'S_m_EG_5': 0.0000, 'S_n_EG_5': 0.0000,
           'S_l_EG_5c': 0.0000, 'S_m_EG_5c': 0.0000, 'S_n_EG_5c': 0.0000,
           'S_l_EG_5d': 0.0000, 'S_m_EG_5d': 0.0000, 'S_n_EG_5d': 0.0000,
           'S_a_3': 0.0000, 'S_b_3': 0.0000, 'S_c_3': 0.0000, 'S_d_3': 0.0000,
           'S_a_EG_1': 0.0000, 'S_a_EG_1q': 0.0000, 'S_a_EG_1v': 0.0000, 'S_a_EG_1u': 0.0000,
           'S_l_EG_55d': 0.0000, 'S_m_EG_55d': 0.0000, 'S_n_EG_55d': 0.0000,
           'S_f_EG_1u5d': 0.0000, 'S_f_EG_1v5': 0.0000, 
           'S_f_EG_1u5': 0.0000, 'S_f_EG_1v5d': 0.0000,
           'S_f_EG_5c1': 0.0000, 'S_f_EG_5c1q': 0.0000,
           'S_a_EG_11q': 0.0000,'S_a_EG_1v1u': 0.0000,
           'S_g_EG_31q': 0.0000,'S_g_EG_31': 0.0000,
           'S_h_EG_5c3': 0.0000,
           'C_11': 0.0000, 'C_12': 0.0000, 'C_44': 0.0000
	},
    'GeSn': {
            'comment': 'Taken from Song et. al. New J. Phys. 21, 073037 (2019)',
           # Lattice parameters
           'lattice_a0': 0, 
           # E_Gamma (eV)
           'EG_1u': 0.0000, 'EG_1q': 0.0000, 'EG_5d': 0.0000, 
           'EG_1v': 0.0000, 'EG_1': 2.2767, 'EG_5': 0.0000, 
           'EG_3': 0.0000, 'EG_5c': 0.0000, 
           # E_Delta (eV)
           'ED_d': 0.0000, 'ED_c': 0.0000, 'ED_so': -4.9535, 
           'ED_dso': 0.0000, 'ED_prime': 0.0000, 'ED_cd_prime': 0.0000,
           # Matrix element (eV.nm)
           'P_P': -0.3346, 'P_Pd': 1.9328, 'P_P2': -1.2452, 'P_P2d': -0.4220, 
           'P_Px': -0.0055, 'P_Pxd': 2.4647, 'P_P3': 0.0089, 'P_P3d': 0.0000, 
           'P_Pu': -0.0441, 'P_Ps': 0.0000, 'P_P_prime': 0.0000, 'P_P2_prime': 0.0000,
           # Strain parameters (Taken from Ref-3)
           'S_l_EG_5': 0.0000, 'S_m_EG_5': 0.0000, 'S_n_EG_5': 0.0000,
           'S_l_EG_5c': 0.0000, 'S_m_EG_5c': 0.0000, 'S_n_EG_5c': 0.0000,
           'S_l_EG_5d': 0.0000, 'S_m_EG_5d': 0.0000, 'S_n_EG_5d': 0.0000,
           'S_a_3': 0.0000, 'S_b_3': 0.0000, 'S_c_3': 0.0000, 'S_d_3': 0.0000,
           'S_a_EG_1': 0.0000, 'S_a_EG_1q': 0.0000, 'S_a_EG_1v': 0.0000, 'S_a_EG_1u': 0.0000,
           'S_l_EG_55d': 0.0000, 'S_m_EG_55d': 0.0000, 'S_n_EG_55d': 0.0000,
           'S_f_EG_1u5d': 0.0000, 'S_f_EG_1v5': 0.0000, 
           'S_f_EG_1u5': 0.0000, 'S_f_EG_1v5d': 0.000,
           'S_f_EG_5c1': 0.0000, 'S_f_EG_5c1q': 0.0000,
           'S_a_EG_11q': 0.0000,'S_a_EG_1v1u': 0.0000,
           'S_g_EG_31q': 0.0000,'S_g_EG_31': 0.0000,
           'S_h_EG_5c3': 0.0000,
           'C_11': 0.0000, 'C_12': 0.0000, 'C_44': 0.0000
	}
    }
