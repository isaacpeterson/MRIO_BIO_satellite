function IUCN_data_params = load_IUCN_params()
    
    IUCN_data_params = struct();
    IUCN_data_params.overwrite_IUCN_data_object = true;
    IUCN_data_params.save_IUCN_tensors = true;
    IUCN_data_params.read_processed_data_from_file = true;
    
    IUCN_data_params.input_data_filepath = '~/Documents/MATLAB/BIO_SATELLITE/IUCN_input_data/';  %%%%%% set to user defined file_path %%%%%%%%
    IUCN_data_params.output_data_filepath = '~/Documents/MATLAB/BIO_SATELLITE/RedList_2016/';
    IUCN_data_params.IUCN_data_type = 'new';
    IUCN_data_params.system_type = 'EORA';
    IUCN_data_params.tensor_type = 'global';
    IUCN_data_params.include_GHG = true;
    IUCN_data_params.read_threat_classification_from_file = false;
    IUCN_data_params.read_IUCN_countries_from_file = false;
    IUCN_data_params.save_processed_IUCN_data = false;
    IUCN_data_params.HSCPC_sector_num = 6357;
    IUCN_data_params.tensor_threat_type = 'threat_group'; %'threat_group' or 'threat_type'
    IUCN_data_params.processed_IUCN_data_filename = [IUCN_data_params.output_data_filepath, IUCN_data_params.system_type '/processed_IUCN_data_', IUCN_data_params.system_type '.mat'];
    
    IUCN_data_params.EORA_concordance_file_prefix = [IUCN_data_params.input_data_filepath 'EORA_threat_concordances/20140807_GlobalMRIO_Conc_IUCN='];  
    IUCN_data_params.allcountriesflag_filename = [IUCN_data_params.input_data_filepath, 'AllCountriesFlag.mat'];
    IUCN_data_params.UN_to_IUCN_codes_filename = [IUCN_data_params.input_data_filepath 'UN_IUCN_codes.txt'];
    IUCN_data_params.EORA_countries_filename = [IUCN_data_params.input_data_filepath 'IUCNcountries.xlsx'];
    IUCN_data_params.EORA_x_filename = [IUCN_data_params.input_data_filepath, 'x_data_NCOUN_187.txt'];
    IUCN_data_params.HSCPC_x_filename = [IUCN_data_params.input_data_filepath 'GlobalRoot.mat'];
    IUCN_data_params.HSCPC_country_codes_filename = [IUCN_data_params.input_data_filepath 'HSCPC_CountryList_updated.txt'];
    IUCN_data_params.new_IUCN_data_threats_filename = [IUCN_data_params.input_data_filepath '2016_All_species_threats.txt'];
    IUCN_data_params.new_IUCN_data_species_filename = [IUCN_data_params.input_data_filepath '2016_All_species_CoO_updated.txt'];
    IUCN_data_params.old_IUCN_data_filename = [IUCN_data_params.input_data_filepath 'fulldatawbirds_usingqrow_ordered.txt'];
    IUCN_data_params.old_threat_cause_class_filename = [IUCN_data_params.input_data_filepath 'ThreatCauseClassification.csv'];
    IUCN_data_params.HSCPC_concordance_filename = [IUCN_data_params.input_data_filepath, '20161026_GlobalMRIO_Conc_Fl_IUCN-ic.csv'];
    IUCN_data_params.EORA_GHG_filename = [IUCN_data_params.input_data_filepath 'GHG_CO2_EORA.txt'];
    IUCN_data_params.tensor_folder = [IUCN_data_params.output_data_filepath IUCN_data_params.system_type, '/IUCN_tensors/'];
end