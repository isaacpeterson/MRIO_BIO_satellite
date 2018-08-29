analyse_MRIO_params = struct();
analyse_MRIO_params.country_of_interest = 'all';
analyse_MRIO_params.threat_direction = 'consumption_based'; %'production_based' or 'consumption_based'
analyse_MRIO_params.sort_type = 'species_num'; % 'threat_num' or 'species_num'
analyse_MRIO_params.aggregate_type = 'MRIO_species_threat_proportion'; % 'MRIO_species_threat_proportion' or 'MRIO_species_threat_path'
analyse_MRIO_params.status_levels_to_use = 'all'; %{'CR', 'EN',  'LC', 'LR_cd', 'LR_lc', 'LR_nt', 'NT', 'VU'};
analyse_MRIO_params.production_col = 2;
analyse_MRIO_params.consumption_col = 1;
analyse_MRIO_params.datapath = '~/Github/MRIO_BIO_SATELLITE/EORA_outputs/';
analyse_MRIO_params.load_MRIO_objects = true;
analyse_MRIO_params.MRIO_x_filename = 'x_data_NCOUN_187.txt';
analyse_MRIO_params.IUCN_data_object_filename = '~/Github/MRIO_BIO_SATELLITE/IUCN_input_data/IUCN_data_object_for_manfred.mat';
analyse_MRIO_params.satellite_species_characteristics_filename = 'satellite_species_characteristics.mat';
analyse_MRIO_params.output_folder = '~/Github/MRIO_BIO_SATELLITE/EORA_outputs/';
analyse_MRIO_params.MRIO_threat_tensor_filename = 'MRIO_threat_tensor.mat';
analyse_MRIO_params.build_threat_tensor = true;
analyse_MRIO_params.rank_type = 'by_industry';

%load(analyse_MRIO_params.IUCN_data_object_filename)  

analyse_MRIO_output_routines(analyse_MRIO_params, IUCN_data_object);