satellite_params = struct();
satellite_params.build_IUCN_data_object = false;
satellite_params.build_IUCN_satellite = false;
satellite_params.display_satellite = false;
satellite_params.collapse_satellite = true;
satellite_params.build_domestic_satellite = false;
satellite_params.build_global_satellite = false;
satellite_params.satellite_collapse_concordance_filename = 'HSCPC_Eora25_secagg.csv';

satellite_params.satellite_collapse_type = 'Eora25';
satellite_params.system_type = 'HSCPC';
satellite_params.IUCN_data_object_filename = [satellite_params.system_type '_IUCN_data_object_2017.mat'];

satellite_params.return_satellite = true;
satellite_params.write_satellite_to_disk = true;

satellite_params.domestic_threats_to_aggregate = 'all';
satellite_params.global_threats_to_aggregate = 'all';

satellite_params.status_levels_to_use = {'CR', 'EN',  'LC', 'DD', 'LR_cd', 'LR_lc', 'LR_nt', 'NT', 'VU'};
satellite_params.country_sort_type = 'EORA';

satellite_params.use_endemics = false;

satellite_params.display_domestic_satellite = true;
satellite_params.display_global_satellite = true;
satellite_params.display_total_satellite = false;
satellite_params.species_sort_type = 'species_class';     %'species_class'
satellite_params.collapse_through_species_sort_type = false;
satellite_params.species_sub_category_type = satellite_params.species_sort_type;
satellite_params.species_sub_category_to_use = {'all'};  %'all' or species classes in capitals
satellite_params.sectors_to_label = 'all'; %{'Colombia'; 'Italy'; 'Finland'; 'Brazil';'Peru';'South Africa';'Madagascar';'Borneo'};
satellite_params.full_species_labels = {'ACTINOPTERYGII','AMPHIBIA','ANTHOZOA',  'ARACHNIDA', 'AVES', 'BIVALVIA', 'CEPHALASPIDOMORPHI', 'CEPHALOPODA', 'CHILOPODA', 'CHONDRICHTHYES', 'CLITELLATA', ...
                                        'DIPLOPODA','ENOPLA','ENTOGNATHA','GASTROPODA','HOLOTHUROIDEA','HYDROZOA','INSECTA','MALACOSTRACA','MAMMALIA', 'MAXILLOPODA','MEROSTOMATA','MYXINI','ONYCHOPHORA','REPTILIA','SARCOPTERYGII'};
satellite_params.species_to_label = satellite_params.full_species_labels; 
satellite_params.overwrite_IUCN_data_object = true;
satellite_params.save_IUCN_tensors = true;
satellite_params.read_processed_data_from_file = true;
satellite_params.input_data_filepath = '~/Github/MRIO_BIO_SATELLITE/IUCN_input_data/';  %%%%%% set to user defined file_path %%%%%%%%
satellite_params.output_data_filepath = '~/Github/MRIO_BIO_SATELLITE/RedList_2016/';
satellite_params.IUCN_data_type = 'new';

satellite_params.tensor_scale = 'country';
satellite_params.include_GHG = true;
satellite_params.read_threat_classification_from_file = false;
satellite_params.read_IUCN_countries_from_file = false;
satellite_params.save_processed_IUCN_data = false;
satellite_params.HSCPC_sector_num = 6357;
satellite_params.tensor_threat_type = 'threat_group'; %'threat_group' or 'threat_type'
satellite_params.satellite_path = [satellite_params.output_data_filepath, satellite_params.system_type, '/', satellite_params.system_type, '_satellite_files/'];
satellite_params.HSCPC_to_Eora_concordance_filepath = '~/Github/MRIO_BIO_SATELLITE/HSCPC_to_Eora_concs/';
satellite_params.EORA_concordance_filepath = [satellite_params.input_data_filepath 'EORA_threat_concordances/'];  
satellite_params.EORA_concordance_file_prefix = '20140807_GlobalMRIO_Conc_IUCN=';  
satellite_params.allcountriesflag_filename = [satellite_params.input_data_filepath, 'AllCountriesFlag.mat'];
satellite_params.UN_to_IUCN_codes_filename = [satellite_params.input_data_filepath 'UN_IUCN_codes.txt'];
satellite_params.EORA_countries_filename = [satellite_params.input_data_filepath 'IUCNcountries.xlsx'];
satellite_params.EORA_x_filename = [satellite_params.input_data_filepath, 'x_data_NCOUN_187.txt'];
satellite_params.HSCPC_x_filename = [satellite_params.input_data_filepath 'GlobalRoot.mat'];
satellite_params.HSCPC_country_codes_filename = [satellite_params.input_data_filepath 'HSCPC_CountryList_updated.txt'];
satellite_params.new_IUCN_data_threats_filename = [satellite_params.input_data_filepath '2016_All_species_threats.txt'];
satellite_params.new_IUCN_data_species_filename = [satellite_params.input_data_filepath '2016_All_species_CoO_updated.txt'];
satellite_params.old_IUCN_data_filename = [satellite_params.input_data_filepath 'fulldatawbirds_usingqrow_ordered.txt'];
satellite_params.old_threat_cause_class_filename = [satellite_params.input_data_filepath 'ThreatCauseClassification.txt'];
satellite_params.HSCPC_concordance_filename = [satellite_params.input_data_filepath, '20161026_GlobalMRIO_Conc_Fl_IUCN-ic.csv'];
satellite_params.EORA_GHG_filename = [satellite_params.input_data_filepath 'GHG_CO2_EORA.txt'];
satellite_params.tensor_folder = [satellite_params.output_data_filepath satellite_params.system_type, '/IUCN_tensors/'];
satellite_params.species_taxons_to_use = 'all';
satellite_params.display_type = 'global';
satellite_params.output_file_type = 'mat';

%IUCN_data_object = process_IUCN_data_routines(satellite_params);
satellite_object = build_IUCN_satellite_routines(IUCN_data_object, satellite_params);