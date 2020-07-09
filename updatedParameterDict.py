HoTCoffeehParameters = {
    'qtnpts'                                    :                                           1,
    'SV_resonanceThreshold'                     :                                         0.0,
    'use_log_fit'                               :                                           1,
    'include_delta_f'                           :                                           1,
    'qynpts'                                    :                                           1,
    'qxnpts'                                    :                                           1,
    'CF_resonanceThreshold'                     :                                         0.0,
    'qznpts'                                    :                                           1,
}

controlParameterList = {
    'niceness'                                  :                                           0,
    'simulation_type'                           :                                    'hybrid',
}

hydroParameters = {
    'Tdec'                                      :                                        0.12,
    'IviscousEqsType'                           :                                           1,
    'Edec'                                      :                                         0.3,
    'ndx'                                       :                                           2,
    'InitialURead'                              :                                           0,
    'IvisBulkFlag'                              :                                           0,
    'dy'                                        :                                         0.1,
    'Initialpitensor'                           :                                           1,
    'T0'                                        :                                         0.6,
    'Ivisflag'                                  :                                           0,
    'ndy'                                       :                                           2,
    'vis'                                       :                                        0.08,
    'IhydroJetoutput'                           :                                           0,
    'iLS'                                       :                                         130,
    'dx'                                        :                                         0.1,
    'visbulknorm'                               :                                         0.0,
    'dt'                                        :                                        0.02,
    'ndt'                                       :                                           5,
}

iSParameters = {
    'turn_on_shear'                             :                                           1,
    'turn_on_bulk'                              :                                           0,
}

iSSParameters = {
    'sample_y_minus_eta_s_range'                :                                         2.0,
    'include_deltaf_bulk'                       :                                           0,
    'number_of_repeated_sampling'               :                                       10000,
    'y_LB'                                      :                                        -2.5,
    'include_deltaf_shear'                      :                                           1,
    'calculate_vn'                              :                                           0,
    'y_RB'                                      :                                         2.5,
    'MC_sampling'                               :                                           2,
    'turn_on_bulk'                              :                                           0,
}

initial_condition_control = {
    'pre-generated_initial_file_path'           :                                'smooth_ICs',
    'pre-generated_initial_file_pattern'        :              'sdAvg_order_[0-9]*_block.dat',
    'pre-generated_initial_file_read_in_mode'   :                                           2,
    'cut_type'                                  :                             'total_entropy',
    'initial_condition_type'                    :                                   'superMC',
    'centrality'                                :                                     '0-10%',
}

photonEmissionParameters = {
    'T_sw_low'                                  :                                      0.1795,
    'T_sw_high'                                 :                                        0.18,
    'enable_polyakov_suppression'               :                                           0,
    'tau_end'                                   :                                        30.0,
    'calHGIdFlag'                               :                                           0,
    'dx'                                        :                                         0.3,
    'dy'                                        :                                         0.3,
    'T_cutlow'                                  :                                         0.1,
    'T_dec'                                     :                                        0.15,
    'differential_flag'                         :                                           0,
    'T_cuthigh'                                 :                                         0.8,
    'dTau'                                      :                                         0.3,
    'tau_start'                                 :                                         0.6,
}

preEquilibriumParameters = {
    'taumin'                                    :                                         0.6,
    'dtau'                                      :                                         0.2,
    'event_mode'                                :                                           1,
    'taumax'                                    :                                         0.6,
}

superMCParameters = {
    'maxx'                                      :                                        13.0,
    'maxy'                                      :                                        13.0,
    'Aproj'                                     :                                         197,
    'ecm'                                       :                                         200,
    'cc_fluctuation_model'                      :                                           6,
    'include_NN_correlation'                    :                                           1,
    'cc_fluctuation_Gamma_theta'                :                                        0.75,
    'finalFactor'                               :                                      28.656,
    'dx'                                        :                                         0.1,
    'dy'                                        :                                         0.1,
    'Atarg'                                     :                                         197,
    'alpha'                                     :                                        0.14,
    'operation'                                 :                                           1,
    'nev'                                       :                                           1,
    'model_name'                                :                                     'MCKLN',
    'lambda'                                    :                                       0.288,
}

