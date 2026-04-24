from .database import material_database
import numpy as np
#import warnings

## ==============================================================================
class _AlloyParams:
    '''
    The functions in this class calculates the parameters for alloy from their
    binary components.
    '''
    _alloy_name_map = {'AlNGaN': 'AlGaN','GaNAlN': 'AlGaN', 'AlGaN':'AlGaN',
                       'SiGe': 'SiGe', 'GeSi': 'SiGe',
                       'SnGe': 'GeSn', 'GeSn': 'GeSn',
                       'SiSn': 'SiSn', 'SnSi': 'SiSn',}
    def __init__(self, binaries=['Si', 'Ge'], alloy_crystal_structure:str='wz', 
                 alloy_type:str=None, use_this_params:dict=None):
        """
        Initiation function of the class _AlloyParams.

        Parameters
        ----------
        binaries : list of strings (case sensitive), optional
            Name of the corresponding binaries of requested alloy. They should
            match the names in database. All implemented materials name list 
            can be found in the README. Alloy 'compositions' will be calculated  
            considering compositions of binaries correspond to x, y, 1-x-y etc. 
            from left to right in the list. For alloy_type='CatAni', it will be counted as
            components it is counted as x, 1-x, y, 1-y etc. 
            The default is ['Si', 'Ge'].
        alloy_crystal_structure :  str, optional [options: 'WZ', 'ZB', 'DM']
            The crystal type of the materials. This will be considered when calculating
            parameters like Poisson ratio etc.
            Use following abbreviation name:
                for wurtzite use 'WZ' or 'wz'.
                for zincblende use 'ZB' or 'zb'.
                for diamond use 'DM' or 'dm'.
            The default is 'wz'. 
        alloy_type : string (case sensitive), optional [options: 'CatAni']
            The alloy type name. Case sensitive. Only needed if alloy is of AxB1-xCxD1-y
            kind. Will be ignored for alloy of type AxB1-x, AxByC1-x-y, AxByCzD1-x-y-z etc.
            The default is None.
        use_this_params : dict
            To use different materials parameters from that given in the database.
            Simply join the binary names to construct the alloy name. e.g.,
            for binaries=['Si', 'Ge'] the alloy name is 'SiGe' or 'GeSi'.
            Material parameters units should be same as in the database.
            e.g. use_mat_params = {'Si': {'P1': 3000}}

        Returns
        -------
        None.

        """
        self.bins_ = list(binaries)
        self.alloy_type_ = alloy_type
        self.alloy_crys_type_ = alloy_crystal_structure.lower()
        # Get parameters from database and update it if use_mat_params is not none.
        self.bin_params_dbs = self._get_params_from_database(use_this_params)
        # Initializing the apply strain keyword
        self.apply_strain_ = False
    
    @staticmethod
    def _update_material_params_locally(use_mat_params, params_db):
        """
        This function allows to update the materials parameters locally before calculation.
        It will not change the original database that comes with the package.

        Parameters
        ----------
        use_mat_params : dict
            To use different materials parameters from that given in the database.
            Simply join the binary names to construct the alloy name. e.g.,
            for binaries=['Si', 'Ge'] the alloy name is 'SiGe' or 'GeSi'.
            Material parameters units should be same as in the database.
            e.g. use_mat_params = {'Si': {'P1': 3000}}
        params_db : dict
            The original material parameters as read from database.

        Returns
        -------
        params_db_tmp : dict
            The updated material parameters after update.

        """
        params_db_tmp = params_db.copy()
        # Avoids appending unnecessary parameters in the return database
        allowed_params = material_database['Si'].keys() # We know one 'default' material in the database
        for mat_name, parms_ in use_mat_params.items():
            if mat_name in params_db_tmp:
                for pms_n, pms_vals in parms_.items():
                    if pms_n in allowed_params:
                        params_db_tmp[mat_name][pms_n] = pms_vals
                    else:
                        print(f"Warning: Ignoring update for {mat_name}:{pms_n}. Does not exist in database.")
            else:
                print(f"Warning: Ignoring update for {mat_name}. Does not exist in loaded database for specific sample.")
                # warnings.warn(f"Ignoring update for {mat_name} material, as it does not exits in the database.", 
                #               RuntimeWarning, stacklevel=1) 
        return params_db_tmp
        
    def _get_params_from_database(self, use_mat_params):
        """
        This function reads requested materials parameters from database and
        updates them if needed.

        Parameters
        ----------
        use_mat_params : dict
            To use different materials parameters from that given in the database.
            Simply join the binary names to construct the alloy name. e.g.,
            for binaries=['Si', 'Ge'] the alloy name is 'SiGe' or 'GeSi'.
            Material parameters units should be same as in the database.
            e.g. use_mat_params = {'Si': {'P1': 3000}}

        Returns
        -------
        dict
            Material parameters.

        """
        alloy_name = ''.join(self.bins_)
        self.alloy_name = _AlloyParams._alloy_name_map.get(alloy_name)
        if self.alloy_name is None: 
            raise ValueError('Requested alloyes are not implemented yet. Contact developer.')
            
        params_db = {self.alloy_name: material_database.get(self.alloy_name).copy()} # These are bowing parameters
        for mat_name in self.bins_:
            params_db[mat_name] = material_database.get(mat_name).copy()

        if use_mat_params is None:
            return params_db  
        else:
            #####_update_database_data_locally
            if alloy_name in use_mat_params:
                use_mat_params[self.alloy_name] = use_mat_params.pop(alloy_name)
            return self._update_material_params_locally(use_mat_params, params_db)
        
    #--------------------------------------------------------------------------
    @staticmethod        
    def _get_substrate_properties(substrate_name):
        """
        Generate the substrate properties for phsedomorphic strain.

        Parameters
        ----------
        substrate_name : str
            The name of the substrate. The name should be in the database. If 
            the name does not exists in the database return None.

        Returns
        -------
        Dictionary
            The parameters of the substrate. Get from database. If substrate
            name does not exists in the database return None.

        """
        return material_database.get(substrate_name).copy()
            
    def _cal_pseudomorphic_strain(self, overwrite_epsilon_in_plane_DCS:float|list|np.ndarray=None):
        a_xs = self.alloy_params_.get('lattice_a0')
        len_axs = len(a_xs)
        
        if overwrite_epsilon_in_plane_DCS is not None:
            if isinstance(overwrite_epsilon_in_plane_DCS, (float, int)):
                epsilon_in_plane_DCS = np.array([overwrite_epsilon_in_plane_DCS]*len_axs)
            else:
                tmp = np.array(overwrite_epsilon_in_plane_DCS)
                if len(tmp) < len_axs:
                    raise ValueError('Overwrite_strain array length should be >= compositions array. Alternatively pass a single float.')
                epsilon_in_plane_DCS = tmp[:len_axs]                
        else:
            if isinstance(self.biaxial_substrate, str):
                substrate_params_dic = self._get_substrate_properties(self.biaxial_substrate)
                substrate_lp = substrate_params_dic.get('lattice_a0') # substrate in-plane lattice parameter
            else:    
                substrate_lp = float(self.biaxial_substrate)
            # Device coordinate system  
            epsilon_in_plane_DCS = substrate_lp/a_xs - 1.0 
        
        if self.alloy_crys_type_ == 'wz':
            #lattice_c = self.alloy_params_.get('lattice_c0') 
            biaxial_dist_coeff = -2*(self.alloy_params_.get('C_13')/self.alloy_params_.get('C_33'))
        elif self.alloy_crys_type_ == 'zb':
            # Crystal coordinate system
            # strain_tensor = [e_xx, e_yy, e_zz, e_yz, e_xz, e_xy]
            # TODO: generalize it using matrix multiplication
            hkl = f'{self.growth_hkl[0]}{self.growth_hkl[1]}{self.growth_hkl[2]}'
            if hkl == '001': # DCS = CCS
                print(f'In-plane strain = {", ".join([f"{x*100:0.2f}%" for x in epsilon_in_plane_DCS])}')
                epsilon_outof_plane_DCS = self._biaxial_Poisson_ratio_zb(hkl)*epsilon_in_plane_DCS
                zeros = np.zeros(np.shape(epsilon_outof_plane_DCS))
                return np.array([epsilon_in_plane_DCS, epsilon_in_plane_DCS, 
                                 epsilon_outof_plane_DCS,zeros,zeros,zeros])
            elif hkl == '100': # DCS = CCS
                epsilon_outof_plane_DCS = self._biaxial_Poisson_ratio_zb(hkl)*epsilon_in_plane_DCS
                zeros = np.zeros(np.shape(epsilon_outof_plane_DCS))
                return np.array([epsilon_outof_plane_DCS, epsilon_in_plane_DCS, 
                                 epsilon_in_plane_DCS,zeros,zeros,zeros])
            elif hkl == '111': # DCS /= CCS
                epsilon_outof_plane_DCS = self._biaxial_Poisson_ratio_zb(hkl)*epsilon_in_plane_DCS
                e_ii = (epsilon_outof_plane_DCS + 2.0*epsilon_in_plane_DCS) / 3.0
                e_ij = (epsilon_outof_plane_DCS - epsilon_in_plane_DCS) / 3.0
                return np.array([e_ii, e_ii, e_ii, e_ij, e_ij, e_ij])
            elif hkl == '110': # DCS /= CCS
                epsilon_outof_plane_DCS = self._biaxial_Poisson_ratio_zb(hkl)*epsilon_in_plane_DCS
                zeros = np.zeros(np.shape(epsilon_outof_plane_DCS))
                e_xxyy = (epsilon_outof_plane_DCS + epsilon_in_plane_DCS) / 2.0
                e_xy = (epsilon_outof_plane_DCS - epsilon_in_plane_DCS) / 2.0
                return np.array([e_xxyy, e_xxyy, epsilon_in_plane_DCS, zeros, zeros, e_xy])
            else:
                raise ValueError('Requested growth direction is not implemented yet. Contact developer.')
        else:
            raise ValueError('Requested alloy crystal structure is not supported yet. Contact developer.')
        return
    
    def _biaxial_Poisson_ratio_zb(self, hkl):
        if hkl in ['001','100','010']:
            # biaxial Poisson ratio = -2(C_12/C_11)
            return -2.0*(self.alloy_params_.get('C_12')/self.alloy_params_.get('C_11'))
        elif hkl == '111': # DCS /= CCS
            # biaxial Poisson ratio = -2[(C_11+2C_12-2C_44)/(C_11+2C_12+4C_44)]
            return -2.0*((self.alloy_params_.get('C_11')+
                          2.0*self.alloy_params_.get('C_12')-
                          2.0*self.alloy_params_.get('C_44'))/
                         (self.alloy_params_.get('C_11')+
                          2.0*self.alloy_params_.get('C_12')+
                          4.0*self.alloy_params_.get('C_44')))
        elif hkl == '110': # DCS /= CCS
            # biaxial Poisson ratio = -[(C_11+3C_12-2C_44)/(C_11+C_12+2C_44)]
            return -1.0*((self.alloy_params_.get('C_11')+
                          3.0*self.alloy_params_.get('C_12')-
                          2.0*self.alloy_params_.get('C_44'))/
                         (self.alloy_params_.get('C_11')+
                          self.alloy_params_.get('C_12')+
                          2.0*self.alloy_params_.get('C_44')))
        
    def _rotation_matrix_from_hkl(hkl:np.ndarray):
        # z-axis surface normal
        z = hkl / np.linalg.norm(hkl)
        return
    #--------------------------------------------------------------------------
    def _get_two_comp_component_alloy_params(self):
        """
        This function calculates the parameters for a 2-composition component alloy 
        from itsbinary component parameters using quadratic interpolation.
        E.g. for any parameter, P:
            P_SixGe1-x = x*P_Si + (1-x)*P_Ge - x*(1-x)*P_bowing 
            P_bowing is the quadratic bowing parameter for the parameter P.
        Returns
        -------
        Parameters for alloy.

        """        
        bin_1_params_db = self.bin_params_dbs.get(self.bins_[0])
        bin_2_params_db = self.bin_params_dbs.get(self.bins_[1])
        alloy_params_db = self.bin_params_dbs.get(self.alloy_name)

        self.alloy_params_ = {}
        for key, bowing in alloy_params_db.items():
            self.alloy_params_[key] = self.comps_ * bin_1_params_db.get(key) +\
                                        (1.0-self.comps_) * bin_2_params_db.get(key)\
                                            - bowing*self.comps_*(1.0-self.comps_)  
        return 
            
    def _get_alloy_params(self, composition:float|np.ndarray):
        """
        This function calculates material parameters for alloy from its
        binary component parameters using interpolation.

        Parameters
        ----------
        

        Returns
        -------
        Parameters for alloy.

        """
        self.comps_ = composition

        if len(self.bins_) == 2:
            self._get_two_comp_component_alloy_params()
        else:
            raise ValueError('Requested multiple composition component alloys are not implemented yet. Contact developer.')
            