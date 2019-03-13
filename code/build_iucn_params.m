function params_object = build_iucn_params()
    params_object = struct();
    params_object.satellite_params = build_satellite_params();
    params_object.build_footprint_params = build_mrio_footprint_params();
    params_object.footprint_analysis_params = build_process_MRIO_params();
end

function mrio_footprint_params = build_mrio_footprint_params()

    mrio_footprint_params = struct();
    mrio_footprint_params.build_mrio_outputs = true;
    mrio_footprint_params.startyear = 2013;
    mrio_footprint_params.endyear = 2013;
    mrio_footprint_params.thresh = 0.01;
    mrio_footprint_params.phase = 199;
    mrio_footprint_params.loop = 82;
    mrio_footprint_params.assessment_scale = 'finalsale';
    mrio_footprint_params.satellite_type = 'global';
    mrio_footprint_params.mrio_output_filename = 'SpThrVals_domestic.mat';
    
    if strcmp(computer,'MACI64');
        mrio_footprint_params.mrio_dir = '~/GitHub/mrio_bio_satellite/processed_data/eora/eora_outputs/';
        mrio_footprint_params.workdir = mrio_footprint_params.mrio_dir;
        mrio_footprint_params.rawdatadir = mrio_footprint_params.mrio_dir;
        mrio_footprint_params.mrio_outputs_file_prefix = [mrio_footprint_params.mrio_dir '21160000_annetest_AllCountries_' sprintf('%03d', mrio_footprint_params.phase)];
        mrio_footprint_params.mrio_outputs_file_suffix = ['_' sprintf('%03d', mrio_footprint_params.loop) '_Markup001.spbin'];
        
    elseif strcmp(computer,'GLNXA64');
    
        mrio_footprint_params.mrio_dir = ['/import/laika1/isa/Eora/Phase' sprintf('%03d', mrio_footprint_params.phase) '/Loop' sprintf('%03d', mrio_footprint_params.loop) '/Results/'];
        mrio_footprint_params.workdir = '/import/emily1/isa/Projects/2018_IUCN_update/';
        mrio_footprint_params.rawdatadir = '/import/emily1/isa/Projects/2017_Tourism/FutuTourism/RawData/'; 
        mrio_footprint_params.concdir = '/import/emily1/isa/Projects/2017_Tourism/FutuTourism/Concordances/';
        mrio_footprint_params.mrio_outputs_file_prefix = [mrio_footprint_params.mrio_dir '21160000_annetest_AllCountries_' sprintf('%03d', mrio_footprint_params.phase)];
        mrio_footprint_params.mrio_outputs_file_suffix = ['_' sprintf('%03d', mrio_footprint_params.loop) '_Markup001.spbin'];
    end
    
    mrio_footprint_params.treedir = [mrio_footprint_params.rawdatadir 'TreeStructure/'];
    mrio_footprint_params.use_sparse_representation = false;
    mrio_footprint_params.satellite_filename = [mrio_footprint_params.workdir 'total_global_satellite.mat'];
    
    mrio_footprint_params.sp_intervals = 1000;
    mrio_footprint_params.q_vec = 1:41000;
    mrio_footprint_params.country_vec = 1:189;
    
    mrio_footprint_params.production_country_exclusions = {'Cayman Islands'};
    mrio_footprint_params.finalsale_country_exclusions = {'Cayman Islands'};
    mrio_footprint_params.production_sector_exclusions = {'Other business activities','Recycling','Business services','Construction of Communication Facilities','Passenger car trade and repairs, service stations','Life Insurnce Service','Toys and games','Other services','Other Services','other services (private)','Technical services for agriculture, forestry, livestock and fishing','Finacial Intermediation and Business Activities'};
    mrio_footprint_params.finalsale_sector_exclusions = {'Agriculture','Services related to agriculture and forestry','Recycling','Crude oil','Life Insurnce Service','Logging','Agriculture, hunting, forestry and fishing'};
    mrio_footprint_params.sector_labels = load('~/GitHub/mrio_bio_satellite/iucn_input_data/EoraFullLabels.mat');

end

function satellite_params = build_satellite_params()

    satellite_params = struct();

    satellite_params.system_type = 'eora';
    satellite_params.build_processed_iucn_data = false;
    satellite_params.overwrite_tensors = true;
    satellite_params.display_satellite = false;
    satellite_params.build_domestic_satellite = true;
    satellite_params.build_global_satellite = true;
    satellite_params.save_processed_iucn_data = true;
    satellite_params.return_satellite = true;
    satellite_params.write_satellite_to_disk = true;
    satellite_params.satellite_collapse_type = 'none';

%     satellite_params.raw_datapath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/';
%     satellite_params.processed_datapath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/';
%     satellite_params.eora_concordance_filepath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/concordances/eora_threat_concordances/';
%     satellite_params.hscpc_to_eora_concordance_filepath = '/import/emily1/isa/IELab/GlobalIELab/RawDataRepository/iucn_redlist/2016/concordances/hscpc_to_eora_concs/';

    satellite_params.raw_datapath = '~/Github/mrio_bio_satellite/iucn_input_data/';
    satellite_params.processed_datapath = '~/Github/mrio_bio_satellite/processed_data/';
    satellite_params.hscpc_to_eora_concordance_filepath = '~/GitHub/mrio_bio_satellite/concordances/hscpc_to_eora_concs/';
    satellite_params.eora_concordance_filepath = '~/GitHub/mrio_bio_satellite/concordances/eora_threat_concordances/';

    satellite_params.satellite_filepath = [satellite_params.processed_datapath, satellite_params.system_type, '/', satellite_params.system_type, '_satellite_files/'];
    satellite_params.pre_processed_tensor_filepath = [satellite_params.processed_datapath satellite_params.system_type, '/', satellite_params.system_type, '_pre_processed_iucn_tensors/'];

    satellite_params.eora_concordance_file_prefix = [satellite_params.eora_concordance_filepath '20140807_GlobalMRIO_Conc_IUCN=']; 
    satellite_params.pre_processed_tensor_file_prefix = [satellite_params.pre_processed_tensor_filepath, satellite_params.system_type, '_iucn_tensor_'];

    satellite_params.processed_iucn_data_filename = [satellite_params.processed_datapath, satellite_params.system_type, '/iucn_data_object_for_manfred.mat'];
    satellite_params.allcountriesflag_filename = [satellite_params.raw_datapath, 'allcountriesflag.mat'];
    satellite_params.un_to_iucn_codes_filename = [satellite_params.raw_datapath 'un_iucn_codes.txt'];
    satellite_params.eora_countries_filename = [satellite_params.raw_datapath 'iucn_countries.xlsx'];
    satellite_params.eora_x_filename = [satellite_params.raw_datapath, 'x_data_187.txt'];
    satellite_params.hscpc_x_filename = [satellite_params.raw_datapath 'GlobalRoot.mat'];
    satellite_params.hscpc_country_codes_filename = [satellite_params.raw_datapath 'hscpc_countrylist_updated.txt'];
    satellite_params.new_iucn_data_threats_filename = [satellite_params.raw_datapath '2016_all_species_threats.txt'];
    satellite_params.new_iucn_data_species_filename = [satellite_params.raw_datapath '2016_all_species_coo_updated.txt'];
    satellite_params.old_iucn_data_filename = [satellite_params.raw_datapath 'fulldatawbirds_usingqrow_ordered.txt'];
    satellite_params.old_threat_cause_class_filename = [satellite_params.raw_datapath 'threatcauseclassification.txt'];
    satellite_params.hscpc_concordance_filename = [satellite_params.eora_concordance_filepath, '20161026_Globalmrio_Conc_Fl_iucn-ic.csv'];
    satellite_params.output_satellite_dir = [satellite_params.processed_datapath, satellite_params.system_type, '/', satellite_params.system_type, '_satellite_files/'];
    satellite_params.eora_ghg_filename = [satellite_params.raw_datapath 'ghg_co2_eora.txt'];
    satellite_params.satellite_collapse_concordance_filename = [satellite_params.hscpc_to_eora_concordance_filepath 'hscpc_eora25_secagg.csv'];

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
end




function footprint_analysis_params = build_process_MRIO_params()
    
    footprint_analysis_params = struct();
    footprint_analysis_params.country_of_interest = 'global';
    footprint_analysis_params.groups_to_count = {'ANIMALIA', 'PLANTAE'};
    footprint_analysis_params.finalsale_scale = 'global'; %'domestic', 'international', 'global'
    footprint_analysis_params.industry_assessment_type = 'finalsale_based'; %'production_based' or 'finalsale_based'
    footprint_analysis_params.sort_data = false;
    footprint_analysis_params.sort_type = 'species_num'; % 'threat_num' or 'species_num'
    footprint_analysis_params.aggregate_type = 'mrio_species_threat_proportion'; % 'mrio_species_threat_proportion' or 'mrio_species_threat_path'
    footprint_analysis_params.status_levels_to_use = 'all'; %{'CR', 'EN',  'LC', 'LR_cd', 'LR_lc', 'LR_nt', 'NT', 'VU'};
    footprint_analysis_params.production_col = 2;
    footprint_analysis_params.consumption_col = 1;

    %analyse_mrio_params.processed_datapath = '/import/emily1/isa/IELab/Roots/GlobalIELab/Roots/GlobalIELab/ProcessedData/iucn_redlist/2016/';
    %analyse_mrio_params.raw_datapath = '/import/emily1/isa/IELab/Roots/GlobalIELab/Roots/GlobalIELab/

    footprint_analysis_params.raw_datapath = '~/Github/mrio_bio_satellite/iucn_input_data/';
    footprint_analysis_params.processed_datapath = '~/Github/mrio_bio_satellite/processed_data/eora/';

    footprint_analysis_params.mrio_x_filename = [footprint_analysis_params.raw_datapath 'x_data_187.txt'];
    footprint_analysis_params.processed_iucn_data_filename = [footprint_analysis_params.processed_datapath, 'iucn_data_object_for_manfred.mat'];
    footprint_analysis_params.satellite_species_characteristics_filename = [footprint_analysis_params.processed_datapath 'eora_satellite_files/satellite_species_characteristics.mat'];
    footprint_analysis_params.output_folder = [footprint_analysis_params.processed_datapath 'consumption_finalsale_production_outputs/'];

    footprint_analysis_params.low_income_countries_filename = [footprint_analysis_params.raw_datapath 'additional_population_data/developing_countries.txt'];
    footprint_analysis_params.finalsale_subs_filename = [footprint_analysis_params.processed_datapath 'eora_outputs/PostExclMarch/SpThrSubs_domestic_final.mat'];
    footprint_analysis_params.finalsale_vals_filename = [footprint_analysis_params.processed_datapath 'eora_outputs/PostExclMarch/SpThrVals_domestic_final.mat'];
    footprint_analysis_params.consumption_level_subs_filename = [footprint_analysis_params.processed_datapath 'eora_outputs/PostExclNov/SpThrSubs_domestic_final.mat'];
    footprint_analysis_params.consumption_level_vals_filename = [footprint_analysis_params.processed_datapath 'eora_outputs/PostExclNov/SpThrVals_domestic_final.mat'];
    footprint_analysis_params.consumption_level_countries_filename = [footprint_analysis_params.processed_datapath 'eora_outputs/PostExclNov/SpThrCnts_domestic_final.mat'];
    
    footprint_analysis_params.load_mrio_objects = true;
    footprint_analysis_params.build_threat_tensor = true;
    footprint_analysis_params.write_expanded_table = true;
    footprint_analysis_params.write_finalsale_country_ranks = true;
    footprint_analysis_params.write_production_country_ranks = true;
    footprint_analysis_params.write_net_country_ranks = false;
    footprint_analysis_params.write_industry_ranks = true;
    footprint_analysis_params.data_threshold = 0;
    footprint_analysis_params.analysis_type = 'by_country';
    footprint_analysis_params.countries_to_exclude = {'Cayman Islands','Netherlands Antilles', 'Former USSR'};
    
end
