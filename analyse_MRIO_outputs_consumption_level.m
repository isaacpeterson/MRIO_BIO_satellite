analyse_MRIO_params = struct();
analyse_MRIO_params.country_of_interest = 'global';
analyse_MRIO_params.assessment_type = 'consumption_based'; %'production_based' or 'consumption_based'
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
analyse_MRIO_params.assessment_scale = 'international'; %'domestic', 'international', 'global'
analyse_MRIO_params.final_sale = false;
analyse_MRIO_params.consumption_level_filename = [analyse_MRIO_params.datapath 'SpThrList_domestic_consumption.mat'] ;

load(analyse_MRIO_params.IUCN_data_object_filename)  

[kingdoms, ~, kingdom_indexes] = unique(RankListShort(:, 1), 'stable');

fid = fopen([analyse_MRIO_params.datapath, analyse_MRIO_params.MRIO_x_filename]);
      x_data = textscan(fid,'%s %s %s %s %f', 'HeaderLines', 0, 'delimiter', ';');
fclose(fid);

unique_countries = unique(x_data{1}, 'stable');
country_indexes = 1:length(unique_countries);
country_index_list = zeros(size(x_data{1}));

for country_ind = country_indexes
	country_index_list(strcmp(x_data{1}, unique_countries(country_ind))) = country_ind;
end
    
industry_characteristics = struct();
industry_characteristics.country_names_list = x_data{1};
industry_characteristics.country_codes_list = x_data{2};
industry_characteristics.country_index_list = country_index_list;
industry_characteristics.country_index_map = [unique_countries num2cell(country_indexes')];
industry_characteristics.commodity_classification_list = x_data{4};

production_country_list = RankListShort(:, 2);
production_commodity_list = RankListShort(:, 3);

final_sale_country_list = RankListShort(:, 4);
final_sale_commodity_list = RankListShort(:, 5);
threat_intensities = RankListShort(indexes_to_use, 6);
consumption_country_list = RankListShort(:, 7);

consumption_country_indexes = cellfun(@(x) find(strcmp(industry_characteristics.country_names_list, x)), consumption_country_list, 'UniformOutput', false);    
production_indexes = cellfun(@(x, y) find(strcmp(industry_characteristics.country_names_list, x) & strcmp(industry_characteristics.commodity_classification_list, y)), production_country_list, production_commodity_list, 'UniformOutput', false);    
final_sale_indexes = cellfun(@(x, y) find(strcmp(industry_characteristics.country_names_list, x) & strcmp(industry_characteristics.commodity_classification_list, y)), final_sale_country_list, final_sale_commodity_list, 'UniformOutput', false);

production_indexes_to_use = cellfun(@(x) length(x) == 1, production_indexes); % quick workaround
final_sale_indexes_to_use = cellfun(@(x) length(x) == 1, final_sale_indexes);
consumption_country_indexes_to_use = cellfun(@(x) length(x) == 1, consumption_country_indexes);
indexes_to_use = production_indexes_to_use & final_sale_indexes_to_use & consumption_country_indexes_to_use;

tensor_indexes = [kingdom_indexes(indexes_to_use) cell2mat(production_indexes(indexes_to_use)) cell2mat(final_sale_indexes(indexes_to_use)), cell2mat(consumption_country_indexes(indexes_to_use))];

aa = sptensor(double(tensor_indexes), double(cell2mat(threat_intensities(indexes_to_use)), double(max(tensor_indexes))));

%analyse_MRIO_output_routines(analyse_MRIO_params, IUCN_data_object);
%name_data = cell2table([IUCN_data_object.UN_to_IUCN_codes.IUCN_country_names IUCN_data_object.UN_to_IUCN_codes.IUCN_industry_codes IUCN_data_object.UN_to_IUCN_codes.UN_industry_codes], 'VariableNames', {'country_name', 'IUCN_code', 'UN_code'});

%writetable(name_data, '~/GitHub/MRIO_BIO_SATELLITE/name_data.txt', 'delimiter', 'tab')
%writetable(t, '~/GitHub/MRIO_BIO_SATELLITE/global_rank_list_with_consumption_country.txt', 'delimiter', 'tab')
 