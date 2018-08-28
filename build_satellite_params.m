function satellite_params = build_satellite_params(IUCN_data_object)

satellite_params = struct();
satellite_params.include_ghg = true;
satellite_params.output_satellite_as_array = true;
satellite_params.write_aggregated_satellite_files = true;
satellite_params.build_domestic_satellite = true;
satellite_params.build_global_satellite = true;
satellite_params.domestic_threats_to_aggregate = 'all';
satellite_params.global_threats_to_aggregate = 'all';
satellite_params.EORA_x_filename = '~/Documents/MATLAB/BIO_SATELLITE/IUCN_input_data/x_data_NCOUN_187.txt';
satellite_params.satellite_filepath = '~/Documents/MATLAB/BIO_SATELLITE/RedList_2016/HSCPC/aggregated_satellite_files/';
satellite_params.IUCN_tensor_filepath = '~/Documents/MATLAB/BIO_SATELLITE/RedList_2016/HSCPC/IUCN_tensors/';
satellite_params.status_levels_to_use = {'CR', 'EN',  'LC', 'DD', 'LR_cd', 'LR_lc', 'LR_nt', 'NT', 'VU'};
satellite_params.country_sort_type = 'EORA';
satellite_params.use_endemics = false;
satellite_params.species_taxons_to_use = IUCN_data_object.IUCN_taxons_list; %IUCN_data_object.IUCN_taxons_list(find(strcmp(IUCN_data_object.country_names_list, 'Australia')));
satellite_params.display_domestic_satellite = true;
satellite_params.display_global_satellite = true;
satellite_params.display_total_satellite = false;
satellite_params.species_sort_type = 'species_class';     %'species_class'
satellite_params.collapse_through_species_sort_type = false;
satellite_params.species_sub_category_type = satellite_params.species_sort_type;
satellite_params.species_sub_category_to_use = {'all'};  %'all' or species classes in capitals
satellite_params.countries_to_label = 'all'; %{'Colombia'; 'Italy'; 'Finland'; 'Brazil';'Peru';'South Africa';'Madagascar';'Borneo'};
satellite_params.full_species_labels = {'ACTINOPTERYGII','AMPHIBIA','ANTHOZOA',  'ARACHNIDA', 'AVES', 'BIVALVIA', 'CEPHALASPIDOMORPHI', 'CEPHALOPODA', 'CHILOPODA', 'CHONDRICHTHYES', 'CLITELLATA', ...
        'DIPLOPODA','ENOPLA','ENTOGNATHA','GASTROPODA','HOLOTHUROIDEA','HYDROZOA','INSECTA','MALACOSTRACA','MAMMALIA', 'MAXILLOPODA','MEROSTOMATA','MYXINI','ONYCHOPHORA','REPTILIA','SARCOPTERYGII'};
satellite_params.species_to_label = satellite_params.full_species_labels; %satellite_params.full_species_labels([1:6, 10, 15:16, 18:21 25:26]);


end