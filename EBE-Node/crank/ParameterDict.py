controlParameterList = {
    'simulation_type'       :   'hydro', 
    # options: 'hybrid', 'hydro', 'hydroEM', 'hydroEM_with_decaycocktail', 
    #          'hydroEM_preEquilibrium', 'hydroEM_with_decaycocktail_with_urqmd'
    'niceness'              :   0,  
    # range from 0 to 19 for process priority, 0 for the highest priority
}

initial_condition_control = {
    'centrality': '0-10%',  # centrality bin
    'cut_type': 'total_entropy',
    # centrality cut variable: total_entropy or Npart
    'initial_condition_type': 'superMC',
    # type of initial conditions: superMC or pre-generated
    'pre-generated_initial_file_path': 'smooth_ICs', 
    # file path for the pre-generated initial condition files
    'pre-generated_initial_file_pattern': 'sdAvg_order_[0-9]*_block.dat',  
    # name pattern for the initial condition files
    'pre-generated_initial_file_read_in_mode': 2, # read in mode for VISH2+1
}

superMCParameters = {
    'model_name'                    :   'MCGlb',    # MCGlb or MCKLN
    'Aproj'                         :   197,
    'Atarg'                         :   197,
    'ecm'                           :   200,
    'finalFactor'                   :   28.656,
    'alpha'                         :   0.14,      # WN/BC mixing ratio in MCGlb
    'lambda'                        :   0.288,      # saturation scale parameter in MCKLN
    'operation'                     :   1,
    'include_NN_correlation'        :   1,
    'cc_fluctuation_model'          :   6,
    'cc_fluctuation_Gamma_theta'    :   0.75,       
    'maxx'                          :   13.0,       # grid size in x (fm)
    'maxy'                          :   13.0,       # grid size in y (fm)
    'dx'                            :   0.1,        # grid spacing in x (fm)
    'dy'                            :   0.1,        # grid spacing in y (fm)
	'nev'							:	1,
}

# only effective when simulation_type == hydroEM_preEquilibrium
preEquilibriumParameters = {
    'event_mode'            :    1,  
    'taumin'                :    0.6,
    'taumax'                :    0.6,
    'dtau'                  :    0.2,
}

hydroParameters = {
    'vis'       	:   0.08,
    'Ivisflag'  :   0,        # flag to use temperature dependent eta/s(T)
    'IvisBulkFlag'  :   1,    # flag for temperature dependence of bulk viscosity
    'visbulknorm'   :   1.0,  # the overall normalization of the bulk viscosity 
                              # (set to 0.0 for shear only simulation)
    'IviscousEqsType'  :  2,  # type of evolution equations for viscous quantities 
                              # (1: Israel-Stewart eq. 2: DNMR eq.)
    'T0'        :   0.6,      # tau_0
    'dt'        :   0.02,     # dtau
    'Tdec'      :   0.12,
    'iLS'       :   130,      # lattice size in transverse plane 2*iLS+1
    'dx'        :   0.10,     # lattice spacing in x (fm) 
                              # need to be the same as dx in superMC
    'dy'        :   0.10,     # lattice spacing in y (fm)
                              # need to be the same as dy in superMC
    'Edec'      :   0.3,  #
    'ndx'       :   2,
    'ndy'       :   2,
    'ndt'       :   5,
    'IhydroJetoutput' :   0,  # switch for output hydro evolution history
    'InitialURead'    :   0,  # set it to be 1 when simulation_type == hydroEM_preEquilibrium
    'Initialpitensor' :   1,  # initialization of pi tensor
                              # 0: initialize to 0; 1: initialze to Navier-Stock value
}

iSSParameters = {
    'turn_on_bulk'                  :   0,
    'include_deltaf_bulk'           :   0,
    'include_deltaf_shear'          :   1,
    'calculate_vn'                  :   0,
    'MC_sampling'                   :   2,
    'number_of_repeated_sampling'   :   10000,
    'y_LB'                          :   -2.5,
    'y_RB'                          :   2.5,
    'sample_y_minus_eta_s_range'    :   2.0,
}

iSParameters = {
    'turn_on_bulk'                  :   0,
    'turn_on_shear'					:   1,
}


photonEmissionParameters = {
    'dx'          :   0.3,
    'dy'          :   0.3,
    'dTau'        :   0.3,
    'tau_start'   :   0.6,
    'tau_end'     :   30.0,
    'T_dec'       :   0.150,
    'T_cuthigh'   :   0.8,
    'T_cutlow'    :   0.1,
    'T_sw_high'   :   0.180,
    'T_sw_low'    :   0.1795,
    'calHGIdFlag' :   0,
    'differential_flag'   :  0,
    'enable_polyakov_suppression'   :    0,
}

HoTCoffeehParameters = {
    'include_delta_f'                   :   1,
    'SV_resonanceThreshold'             :   0.00,
    'CF_resonanceThreshold'             :   0.00,
    'use_log_fit'                       :   1,
    'qtnpts'                            :   1,
    'qxnpts'                            :   1,
    'qynpts'                            :   1,
    'qznpts'                            :   1,
}
