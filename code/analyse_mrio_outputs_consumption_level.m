analyse_mrio_params = struct();
analyse_mrio_params.country_of_interest = 'global';
analyse_mrio_params.groups_to_count = {'ANIMALIA', 'PLANTAE'};
analyse_mrio_params.finalsale_scale = 'global'; %'domestic', 'international', 'global'
analyse_mrio_params.industry_assessment_type = 'finalsale_based'; %'production_based' or 'finalsale_based'
analyse_mrio_params.sort_data = false;
analyse_mrio_params.sort_type = 'species_num'; % 'threat_num' or 'species_num'
analyse_mrio_params.aggregate_type = 'mrio_species_threat_proportion'; % 'mrio_species_threat_proportion' or 'mrio_species_threat_path'
analyse_mrio_params.status_levels_to_use = 'all'; %{'CR', 'EN',  'LC', 'LR_cd', 'LR_lc', 'LR_nt', 'NT', 'VU'};
analyse_mrio_params.production_col = 2;
analyse_mrio_params.consumption_col = 1;
analyse_mrio_params.
analyse_mrio_params.datapath = '~/Github/mrio_bio_satellite/eora_outputs/';
analyse_mrio_params.load_mrio_objects = true;
analyse_mrio_params.mrio_x_filename = 'x_data_187.txt';
analyse_mrio_params.iucn_data_object_filename = '~/Github/mrio_bio_satellite/iucn_input_data/iucn_data_object_for_manfred.mat';
analyse_mrio_params.satellite_species_characteristics_filename = 'satellite_species_characteristics.mat';
analyse_mrio_params.output_folder = '~/Github/mrio_bio_satellite/mrio_output_tables/';
analyse_mrio_params.mrio_threat_tensor_filename = 'mrio_threat_tensor.mat';
analyse_mrio_params.low_income_countries_filename = '~/Github/mrio_bio_satellite/additional_population_data/developing_countries.txt';
analyse_mrio_params.build_threat_tensor = true;
analyse_mrio_params.write_expanded_table = true;
analyse_mrio_params.write_finalsale_country_ranks = true;
analyse_mrio_params.write_production_country_ranks = true;
analyse_mrio_params.write_net_country_ranks = false;
analyse_mrio_params.write_industry_ranks = true;
analyse_mrio_params.finalsale = false;
analyse_mrio_params.data_threshold = 0;
analyse_mrio_params.analysis_type = 'by_country';
analyse_mrio_params.countries_to_exclude = {'Cayman Islands','Netherlands Antilles', 'Former USSR'};

load(analyse_mrio_params.iucn_data_object_filename)  

load([analyse_mrio_params.datapath analyse_mrio_params.satellite_species_characteristics_filename]);
[~, ~, inds_to_use] = intersect(species_characteristics.species_taxons, iucn_data_object.iucn_threat_taxons, 'stable');
species_characteristics.species_kingdom = iucn_data_object.iucn_species_kingdom(inds_to_use);

trade_characteristics = analyse_global_consumption_routines(iucn_data_object, analyse_mrio_params, species_characteristics); 

writetable(trade_characteristics.aggregated_sector_scale.consumption_to_finalsale_table, '~/Github/mrio_bio_satellite/mrio_output_tables/aggregated_sector_scale_consumption_finalsale_table.txt', 'delimiter', 'tab')
writetable(trade_characteristics.country_scale.net_trade_characteristics.table, '~/Github/mrio_bio_satellite/mrio_output_tables/country_scale_net_table.txt', 'delimiter', 'tab')

% find(strcmp(iucn_data_object.industry_characteristics.country_names_list(trade_characteristics.finalsale_data.sector_to_sector_scale.aggregated_paths(:, 1)), 'China')...
% & strcmp(iucn_data_object.industry_characteristics.commodity_classification_list(trade_characteristics.finalsale_data.sector_to_sector_scale.aggregated_paths(:, 1)), 'Construction')...
% & strcmp(iucn_data_object.industry_characteristics.country_names_list(trade_characteristics.finalsale_data.sector_to_sector_scale.aggregated_paths(:, 2)), 'China')...
% & strcmp(iucn_data_object.industry_characteristics.commodity_classification_list(trade_characteristics.finalsale_data.sector_to_sector_scale.aggregated_paths(:, 2)), 'Forestry'))
% 
% china_forestry_species_table = cell2table([species_characteristics.species_names(trade_characteristics.finalsale_data.sector_to_sector_scale.grouped_aggregates{754}) ...
%                                             num2cell(trade_characteristics.finalsale_data.sector_to_sector_scale.grouped_path_vals{754})], ...
%                                             'VariableNames', {'threatened_species', 'threat_intensity'});
% 
% writetable(china_forestry_species_table, '~/Github/mrio_bio_satellite/mrio_output_tables/china_construction_china_forestry_species_list.txt', 'delimiter', 'tab')  

international_indexes = find(~strcmp(trade_characteristics.aggregated_sector_scale.consumption_to_finalsale_table.consumption_country, ...
                                     trade_characteristics.aggregated_sector_scale.consumption_to_finalsale_table.finalsale_country));

T = trade_characteristics.aggregated_sector_scale.consumption_to_finalsale_table([1:10 international_indexes(1:10)'], 1:7);
T.attributed_threat_intensity = round(T.attributed_threat_intensity, 1);
writetable(T, '~/Github/mrio_bio_satellite/mrio_output_tables/aggregated_sector_scale_consumption_finalsale_table_consise.txt', 'delimiter', 'tab')

    