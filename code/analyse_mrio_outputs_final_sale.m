analyse_MRIO_params = struct();
analyse_MRIO_params.country_of_interest = 'global';
analyse_MRIO_params.assessment_scale = 'global'; %'domestic', 'international', 'global'
analyse_MRIO_params.industry_assessment_type = 'final_sale_based'; %'production_based' or 'final_sale_based'
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
analyse_MRIO_params.output_folder = '~/Documents/';
analyse_MRIO_params.MRIO_threat_tensor_filename = 'MRIO_threat_tensor.mat';
analyse_MRIO_params.build_threat_tensor = true;
analyse_MRIO_params.write_expanded_table = true;
analyse_MRIO_params.write_consumption_country_ranks = true;
analyse_MRIO_params.write_production_country_ranks = true;
analyse_MRIO_params.write_net_country_ranks = false;
analyse_MRIO_params.write_industry_ranks = true;
analyse_MRIO_params.final_sale = false;

load([analyse_MRIO_params.datapath '/PostExclMarch/SpThrSubs_domestic_final.mat'] )
load([analyse_MRIO_params.datapath '/PostExclMarch/SpThrVals_domestic_final.mat'] )
load(analyse_MRIO_params.IUCN_data_object_filename) 

MRIO_threat_tensor = sptensor(double(SpThrSubs), double(SpThrVals), double(max(SpThrSubs)));
    
load([analyse_MRIO_params.datapath analyse_MRIO_params.satellite_species_characteristics_filename]);
[~, ~, inds_to_use] = intersect(species_characteristics.species_taxons, IUCN_data_object.IUCN_threat_taxons, 'stable');
species_characteristics.species_kingdom = IUCN_data_object.IUCN_species_kingdom(inds_to_use);
    
final_sale_data = analyse_MRIO_output_routines(analyse_MRIO_params, IUCN_data_object, MRIO_threat_tensor, species_characteristics);

% G12 = {'United States', 'Australia', 'Belgium', 'Canada', 'France', 'Germany', 'Italy', 'Japan', 'Netherlands', 'Spain','Sweden', 'Switzerland', 'United Kingdom'}
 %name_data = cell2table([IUCN_data_object.UN_to_IUCN_codes.IUCN_country_names IUCN_data_object.UN_to_IUCN_codes.IUCN_industry_codes IUCN_data_object.UN_to_IUCN_codes.UN_industry_codes], 'VariableNames', {'country_name', 'IUCN_code', 'UN_code'});

%writetable(name_data, '~/GitHub/MRIO_BIO_SATELLITE/name_data.txt', 'delimiter', 'tab')
%writetable(t, '~/GitHub/MRIO_BIO_SATELLITE/global_rank_list_with_consumption_country.txt', 'delimiter', 'tab')