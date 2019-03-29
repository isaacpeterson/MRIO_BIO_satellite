function params_object = build_iucn_params()

        params_object.global_params = build_global_params;
        params_object.footprint_objects_params = build_footprint_objects_params(params_object.global_params);
        params_object.satellite_params = build_satellite_params(params_object.global_params);
        params_object.build_footprint_params = build_footprint_params(params_object.global_params);
        params_object.analyse_footprint_params = analyse_footprint_params(params_object.global_params);
    
end


function global_params = build_global_params()
    
    global_params.satellite_type = 'direct';
    global_params.system_type = 'eora';
    global_params.phase = 199;
    global_params.loop = 82;

    if strcmp(computer,'MACI64');
        
        global_params.processed_data_filepath = ['~/GitHub/mrio_bio_satellite/processed_data/', global_params.system_type, '/'];
        global_params.mrio_filepath = [global_params.processed_data_filepath 'mrio_outputs/'];
        global_params.processed_mrio_filepath = [global_params.processed_data_filepath 'mrio_outputs/'];
        
        global_params.raw_filepath = '~/Github/mrio_bio_satellite/iucn_input_data/';
        global_params.concordance_filepath = '~/GitHub/mrio_bio_satellite/concordances/eora_threat_concordances/';
        
    elseif strcmp(computer,'GLNXA64');
        
        global_params.processed_data_filepath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/';
        global_params.mrio_filepath = ['/import/laika1/isa/Eora/Phase' sprintf('%03d', global_params.phase) '/Loop' sprintf('%03d', global_params.loop) '/Results/'];
        global_params.processed_mrio_filepath = [global_params.processed_data_filepath 'mrio_outputs/'];
        
        global_params.raw_filepath = '/import/emily1/isa/Projects/2017_Tourism/FutuTourism/RawData/'; 
%         global_params.raw_filepath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/';

        global_params.concordance_filepath = '/import/emily1/isa/Projects/2017_Tourism/FutuTourism/Concordances/';
        
        analyse_mrio_params.processed_filepath = '/import/emily1/isa/IELab/Roots/GlobalIELab/Roots/GlobalIELab/ProcessedData/iucn_redlist/2016/';
        %     global_params.raw_filepath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/';
%     global_params.processed_data_filepath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/';
%     global_params.eora_concordance_filepath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/concordances/';
%     global_params.hscpc_to_eora_concordance_filepath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/concordances/hscpc_to_eora_concs/';
    end
end

function build_footprint_params = build_footprint_params(global_params)

    build_footprint_params = struct();
    build_footprint_params.build_footprint = false;
    build_footprint_params.build_mrio_outputs = false;
    build_footprint_params.thresh = 1e-4;
    
    build_footprint_params.mrio_data_file_prefix = [global_params.mrio_filepath '21160000_annetest_AllCountries_' sprintf('%03d', global_params.phase)];
    build_footprint_params.mrio_data_file_suffix = ['_' sprintf('%03d', global_params.loop) '_Markup001.spbin'];
    build_footprint_params.mrio_objects_file_prefix = [global_params.processed_mrio_filepath 'mrio_objects_'];
    
    build_footprint_params.use_sparse_representation = false;
    build_footprint_params.satellite_filename = [global_params.processed_data_filepath 'collapsed_satellites/' global_params.satellite_type '_threats_satellite.mat'];
    build_footprint_params.countries_to_assess = 'all';
    build_footprint_params.species_to_assess = 'all';
    build_footprint_params.production_country_exclusions = {'Cayman Islands'};
    build_footprint_params.finalsale_country_exclusions = {'Cayman Islands'};
    build_footprint_params.production_sector_exclusions = {'Other business activities','Recycling','Business services','Construction of Communication Facilities','Passenger car trade and repairs, service stations','Life Insurnce Service','Toys and games','Other services','Other Services','other services (private)','Technical services for agriculture, forestry, livestock and fishing','Finacial Intermediation and Business Activities'};
    build_footprint_params.finalsale_sector_exclusions = {'Agriculture','Services related to agriculture and forestry','Recycling','Crude oil','Life Insurnce Service','Logging','Agriculture, hunting, forestry and fishing'};
    build_footprint_params.sector_labels = load('~/GitHub/mrio_bio_satellite/iucn_input_data/EoraFullLabels.mat');
    
    build_footprint_params.footprint_filename_prefix = [global_params.processed_data_filepath 'footprints/'];
    build_footprint_params.footprint_filename_suffix = ['_' global_params.satellite_type, '_footprints_marine.mat'];

end


function satellite_params = build_satellite_params(global_params)

    satellite_params = struct();

    satellite_params.build_processed_iucn_data = true;
    satellite_params.overwrite_tensors = true;
    satellite_params.display_satellite = false;
    satellite_params.save_processed_data = false;
    satellite_params.return_satellite = true;
    satellite_params.return_satellite_array = true;
    satellite_params.write_satellite_to_disk = true;
    satellite_params.satellite_collapse_type = 'none';
    satellite_params.build_direct_satellite = true;
    satellite_params.build_greenhouse_satellite = true;
    satellite_params.satellite_filepath = [global_params.processed_data_filepath, global_params.system_type, '_satellite_files/'];
    satellite_params.eora_concordance_file_prefix = [global_params.concordance_filepath '20140807_GlobalMRIO_Conc_IUCN=']; 
    satellite_params.iucn_tensor_file_prefix = [global_params.processed_data_filepath, global_params.system_type,  '_iucn_tensors/'];

    satellite_params.processed_iucn_data_filename_prefix = [global_params.processed_data_filepath, 'iucn_data_object_'];
    
    satellite_params.greenhouse_flag_filename = [global_params.raw_filepath, 'allcountriesflag.mat'];
    satellite_params.iucn_data_threats_filename = [global_params.raw_filepath '2016_all_species_threats.txt'];
    satellite_params.iucn_data_species_filename = [global_params.raw_filepath '2016_all_species_coo_updated.txt'];
    satellite_params.old_iucn_data_filename = [global_params.raw_filepath 'fulldatawbirds_usingqrow_ordered.txt'];
    satellite_params.old_threat_cause_class_filename = [global_params.raw_filepath 'threatcauseclassification.txt'];
    satellite_params.hscpc_concordance_filename = [global_params.concordance_filepath, '20161026_Globalmrio_Conc_Fl_iucn-ic.csv'];
    satellite_params.output_satellite_dir = [global_params.processed_data_filepath, global_params.system_type, '_satellite_files/'];
    satellite_params.satellite_collapse_concordance_filename = [global_params.concordance_filepath 'hscpc_eora25_secagg.csv'];
    
        
    satellite_params.species_taxons_to_use = 'all';
    satellite_params.display_type = 'global';
    satellite_params.output_file_type = 'mat';
    satellite_params.overwrite_processed_iucn_data = true;
    satellite_params.read_processed_data_from_file = true;
    satellite_params.iucn_data_type = 'new';
    satellite_params.tensor_scale = 'country';
    satellite_params.include_ghg = true;
    satellite_params.read_threat_classification_from_file = false;
    satellite_params.read_iucn_countries_from_file = false;
    satellite_params.hscpc_sector_num = 6357;
    satellite_params.tensor_threat_type = 'threat_group'; %'threat_group' or 'threat_type'
    satellite_params.direct_threats = 'all';
    satellite_params.greenhouse_threats = 'all';
    satellite_params.status_levels_to_use = {'CR', 'EN',  'LC', 'DD', 'LR_cd', 'LR_lc', 'LR_nt', 'NT', 'VU'};
    satellite_params.country_sort_type = global_params.system_type;
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
    
end

function footprint_objects_params = build_footprint_objects_params(global_params)
    footprint_objects_params.satellite_species_characteristics_filename = [global_params.processed_data_filepath, '/satellite_species_characteristics.mat'];
    footprint_objects_params.un_to_iucn_codes_filename = [global_params.raw_filepath 'un_iucn_codes.txt'];
    footprint_objects_params.eora_countries_filename = [global_params.raw_filepath 'iucn_countries.xlsx'];
    footprint_objects_params.eora_x_filename = [global_params.raw_filepath, 'x_data_187.txt'];
    footprint_objects_params.hscpc_x_filename = [global_params.raw_filepath 'GlobalRoot.mat'];
    footprint_objects_params.low_income_countries_filename = [global_params.raw_filepath 'additional_population_data/developing_countries.txt'];
%     footprint_objects_params.industry_labels_filename = [global_params.raw_filepath 'EoraFullLabels.mat'];
    footprint_objects_params.eora_countries_filename = [global_params.raw_filepath 'iucn_countries.xlsx'];
    footprint_objects_params.hscpc_country_codes_filename = [global_params.raw_filepath 'hscpc_countrylist_updated.txt'];
    footprint_objects_params.eora_ghg_filename = [global_params.raw_filepath 'ghg_co2_eora.txt'];
end

function analyse_footprint_params = analyse_footprint_params(global_params)
    
    analyse_footprint_params = struct();
    analyse_footprint_params.scale_ecology_groups = true;
    
    analyse_footprint_params.country_of_interest = 'Australia';
    analyse_footprint_params.sort_data = false;
    analyse_footprint_params.sort_type = 'species_num'; % 'threat_num' or 'species_num'
    analyse_footprint_params.output_folder = [global_params.processed_data_filepath 'consumption_finalsale_production_outputs/'];
    analyse_footprint_params.species_group_weights = [1 1 1 1];
    analyse_footprint_params.weight_species = false;

    analyse_footprint_params.write_expanded_table = true;
    analyse_footprint_params.write_finalsale_country_ranks = true;
    analyse_footprint_params.write_production_country_ranks = true;
    analyse_footprint_params.write_net_country_ranks = false;
    analyse_footprint_params.write_industry_ranks = true;
    analyse_footprint_params.data_threshold = 0;
    analyse_footprint_params.countries_to_exclude = {'Cayman Islands','Netherlands Antilles', 'Former USSR'};
    
end
