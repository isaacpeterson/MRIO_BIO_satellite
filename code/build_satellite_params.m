function satellite_params = build_satellite_params()
satellite_params = struct();
satellite_params.system_type = 'eora';
satellite_params.build_iucn_data_object = false;

satellite_params.display_satellite = false;
satellite_params.build_domestic_satellite = true;
satellite_params.build_global_satellite = true;

satellite_params.save_processed_iucn_data = true;
satellite_params.overwrite_tensors = false;
satellite_params.species_taxons_to_use = 'all';
satellite_params.display_type = 'global';
satellite_params.output_file_type = 'mat';
satellite_params.overwrite_iucn_data_object = true;
satellite_params.read_processed_data_from_file = true;
satellite_params.iucn_data_type = 'new';
satellite_params.tensor_scale = 'country';
satellite_params.include_ghg = true;
satellite_params.read_threat_classification_from_file = false;
satellite_params.read_iucn_countries_from_file = false;
satellite_params.hscpc_sector_num = 6357;
satellite_params.tensor_threat_type = 'threat_group'; %'threat_group' or 'threat_type'
satellite_params.return_satellite = true;
satellite_params.write_satellite_to_disk = true;

satellite_params.satellite_collapse_type = 'none';

satellite_params.domestic_threats_to_aggregate = 'all';
satellite_params.global_threats_to_aggregate = 'all';

satellite_params.status_levels_to_use = {'CR', 'EN',  'LC', 'DD', 'LR_cd', 'LR_lc', 'LR_nt', 'NT', 'VU'};
satellite_params.country_sort_type = 'eora';

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


satellite_params.input_data_filepath = '~/GitHub/mrio_bio_satellite/iucn_input_data/';  %%%%%% set to user defined file_path %%%%%%%%
satellite_params.output_data_filepath = '~/GitHub/mrio_bio_satellite/satellite_accounts/2016/';
satellite_params.satellite_path = [satellite_params.output_data_filepath, satellite_params.system_type, '/', satellite_params.system_type, '_satellite_files/'];
satellite_params.hscpc_to_eora_concordance_filepath = '~/GitHub/mrio_bio_satellite/hscpc_to_eora_concs/';
satellite_params.eora_concordance_file_prefix = [satellite_params.input_data_filepath 'eora_threat_concordances/' '20140807_Globalmrio_Conc_IUCN='];  
satellite_params.allcountriesflag_filename = [satellite_params.input_data_filepath, 'allcountriesflag.mat'];
satellite_params.un_to_iucn_codes_filename = [satellite_params.input_data_filepath 'un_iucn_codes.txt'];
satellite_params.eora_countries_filename = [satellite_params.input_data_filepath 'iucn_countries.xlsx'];
satellite_params.eora_x_filename = [satellite_params.input_data_filepath, 'x_data_187.txt'];
satellite_params.hscpc_x_filename = [satellite_params.input_data_filepath 'GlobalRoot.mat'];
satellite_params.hscpc_country_codes_filename = [satellite_params.input_data_filepath 'hscpc_countrylist_updated.txt'];
satellite_params.new_iucn_data_threats_filename = [satellite_params.input_data_filepath '2016_all_species_threats.txt'];
satellite_params.new_iucn_data_species_filename = [satellite_params.input_data_filepath '2016_all_species_coo_updated.txt'];
satellite_params.old_iucn_data_filename = [satellite_params.input_data_filepath 'fulldatawbirds_usingqrow_ordered.txt'];
satellite_params.old_threat_cause_class_filename = [satellite_params.input_data_filepath 'threatcauseclassification.txt'];
satellite_params.hscpc_concordance_filename = [satellite_params.input_data_filepath, '20161026_Globalmrio_Conc_Fl_iucn-ic.csv'];
satellite_params.output_satellite_dir = [satellite_params.output_data_filepath, satellite_params.system_type, '/', satellite_params.system_type, '_satellite_files/'];
satellite_params.eora_ghg_filename = [satellite_params.input_data_filepath 'ghg_co2_eora.txt'];
satellite_params.tensor_folder = [satellite_params.output_data_filepath satellite_params.system_type, '/', satellite_params.system_type, '_iucn_tensors/'];
satellite_params.iucn_data_object_filename = [satellite_params.system_type '_iucn_data_object_2017.mat'];
satellite_params.satellite_collapse_concordance_filename = [satellite_params.hscpc_to_eora_concordance_filepath 'hscpc_eora25_secagg.csv'];
end