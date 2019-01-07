analyse_MRIO_params = struct();
analyse_MRIO_params.country_of_interest = 'global';
analyse_MRIO_params.groups_to_count = {'ANIMALIA', 'PLANTAE'};
analyse_MRIO_params.assessment_scale = 'global'; %'domestic', 'international', 'global'
analyse_MRIO_params.industry_assessment_type = 'finalsale_based'; %'production_based' or 'finalsale_based'
analyse_MRIO_params.sort_data = false;
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
analyse_MRIO_params.write_finalsale_country_ranks = true;
analyse_MRIO_params.write_production_country_ranks = true;
analyse_MRIO_params.write_net_country_ranks = false;
analyse_MRIO_params.write_industry_ranks = true;
analyse_MRIO_params.finalsale = false;

analyse_MRIO_params.analysis_type = 'by_country';

load(analyse_MRIO_params.IUCN_data_object_filename)  

load([analyse_MRIO_params.datapath analyse_MRIO_params.satellite_species_characteristics_filename]);
[~, ~, inds_to_use] = intersect(species_characteristics.species_taxons, IUCN_data_object.IUCN_threat_taxons, 'stable');
species_characteristics.species_kingdom = IUCN_data_object.IUCN_species_kingdom(inds_to_use);

unique_countries = unique(IUCN_data_object.industry_characteristics.country_names_list, 'stable');

finalsale_subs = load([analyse_MRIO_params.datapath '/PostExclMarch/SpThrSubs_domestic_final.mat'] );
finalsale_vals = load([analyse_MRIO_params.datapath '/PostExclMarch/SpThrVals_domestic_final.mat'] );
finalsale_threat_tensor = sptensor(double(finalsale_subs.SpThrSubs), double(finalsale_vals.SpThrVals), double(max(finalsale_subs.SpThrSubs)));

disp('processing outputs at final sale level')
finalsale_data = analyse_MRIO_output_routines(analyse_MRIO_params, IUCN_data_object, finalsale_threat_tensor, species_characteristics);


%%%%%%%%%%%%%%%% GLOBAL LEVEL ANALYSIS
% global_MRIO_threat_tensor = sptensor(double(SpThrSubs), double(SpThrVals), double(max(SpThrSubs)));
% global_cons_MRIO = analyse_MRIO_output_routines(analyse_MRIO_params, IUCN_data_object, global_MRIO_threat_tensor, species_characteristics);
% 
% [~, sorted_inds] = sort(cell2mat(global_cons_MRIO.ranked_countries_finalsale(:, 4)), 'descend');
% T = cell2table(global_cons_MRIO.ranked_countries_finalsale(sorted_inds, [1 4]), 'VariableNames', {'consumption_country', 'threat_intensity'});
% writetable(T, '~/GitHub/MRIO_BIO_SATELLITE/global_consumption_table_simple.txt', 'delimiter', 'tab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

consumption_subs = load([analyse_MRIO_params.datapath '/PostExclNov/SpThrSubs_domestic_final.mat'] );
consumption_vals = load([analyse_MRIO_params.datapath '/PostExclNov/SpThrVals_domestic_final.mat']);
consumption_countries = load([analyse_MRIO_params.datapath '/PostExclNov/SpThrCnts_domestic_final.mat']);

country_num = max(consumption_countries.SpThrCountries);
consumption_country_set = cell(country_num, 1);

for country_index = 1:country_num
    
    if strcmp(analyse_MRIO_params.analysis_type, 'by_country')    
        current_country_set = (consumption_countries.SpThrCountries == country_index);
    else
        current_country_set = 1:length(consumption_countries.SpThrCountries);
    end

    MRIO_threat_tensor = sptensor(double(consumption_subs.SpThrSubs(current_country_set, :)), double(consumption_vals.SpThrVals(current_country_set)), double(max(consumption_subs.SpThrSubs(current_country_set, :))));
    
    consumption_country_set{country_index} = analyse_MRIO_output_routines(analyse_MRIO_params, IUCN_data_object, MRIO_threat_tensor, species_characteristics);
   
    disp(['country ', num2str(country_index), ' of ', num2str(country_num), ' done'])
    
end
    

global_aggregated_consumption = struct();

%%%%%%%   INDUSTRY LEVEL ROUTINES 

global_aggregated_consumption_per_industry = cellfun(@(x, y) [repmat({y}, [size(x.industry_level.aggregated_path_vals, 1) 1]) ...
                                             IUCN_data_object.industry_characteristics.country_names_list(x.industry_level.aggregated_paths)  ...
                                             IUCN_data_object.industry_characteristics.commodity_classification_list(x.industry_level.aggregated_paths) ...
                                             num2cell(x.industry_level.aggregated_path_vals)], ...
                                             consumption_country_set, unique_countries, 'UniformOutput', false);

global_aggregated_consumption_per_industry = vertcat(global_aggregated_consumption_per_industry{:});

[~, sorted_path_indexes] = sort(cell2mat(global_aggregated_consumption_per_industry(:, end)), 'descend');
global_aggregated_consumption_per_industry_sorted = global_aggregated_consumption_per_industry(sorted_path_indexes, :);

global_aggregated_consumption_industry_index_block = cellfun(@(x, y) ...
[repmat({y}, [size(x.industry_level.aggregated_paths, 1) 1]) num2cell(x.industry_level.aggregated_paths) num2cell(x.industry_level.aggregated_path_vals)], ...
    consumption_country_set, unique_countries, 'UniformOutput', false);


%%%%%%%   INDUSTRY LEVEL ROUTINES 

global_aggregated_consumption.industry_level.country_index_list = cellfun(@(x, y) repmat({y}, [size(x.industry_level.aggregated_paths, 1) 1]), consumption_country_set, num2cell((1:length(unique_countries))'), 'UniformOutput', false);
global_aggregated_consumption.industry_level.country_index_list = cell2mat(vertcat(global_aggregated_consumption.industry_level.country_index_list{:}));

global_aggregated_consumption.industry_level.country_list = cellfun(@(x, y) repmat({y}, [size(x.industry_level.aggregated_paths, 1) 1]), consumption_country_set, unique_countries, 'UniformOutput', false);
global_aggregated_consumption.industry_level.country_list = vertcat(global_aggregated_consumption.industry_level.country_list{:});

global_aggregated_consumption.industry_level.index_list = cellfun(@(x) num2cell(x.industry_level.aggregated_paths), consumption_country_set, 'UniformOutput', false);
global_aggregated_consumption.industry_level.index_list = cell2mat(vertcat(global_aggregated_consumption.industry_level.index_list{:}));

global_aggregated_consumption.industry_level.threat_intensities = cellfun(@(x) num2cell(x.industry_level.aggregated_path_vals), consumption_country_set, 'UniformOutput', false);
global_aggregated_consumption.industry_level.threat_intensities = cell2mat(vertcat(global_aggregated_consumption.industry_level.threat_intensities{:}));


industry_link_indexes = arrayfun(@(x) find(global_aggregated_consumption.industry_level.index_list == x), finalsale_data.industry_level.aggregated_paths, 'UniformOutput', false);

non_empties = cellfun('length', industry_link_indexes) > 0;
industry_link_indexes = industry_link_indexes(non_empties);

grouped_industry_index_links = cellfun(@(x) [global_aggregated_consumption.industry_level.country_list(x)...
                                IUCN_data_object.industry_characteristics.country_names_list(global_aggregated_consumption.industry_level.index_list(x))  ...
                               IUCN_data_object.industry_characteristics.commodity_classification_list(global_aggregated_consumption.industry_level.index_list(x))], industry_link_indexes, 'UniformOutput', false);


scaled_global_industry_threat_intensities = cellfun(@(x, y) [global_aggregated_consumption.industry_level.threat_intensities(x) ./ (sum(global_aggregated_consumption.industry_level.threat_intensities(x))) ...
                                                            global_aggregated_consumption.industry_level.threat_intensities(x) .* (y ./ sum(global_aggregated_consumption.industry_level.threat_intensities(x)))], ...
                                                          industry_link_indexes, num2cell(finalsale_data.industry_level.aggregated_path_vals(non_empties)), 'UniformOutput', false);

           
%finalsale_data_to_use = [finalsale_data.industry_level.aggregates.species_counts{:} finalsale_data.industry_level.aggregated_path_vals];

finalsale_data_to_use = [finalsale_data.industry_level.aggregates.species_counts{:} finalsale_data.industry_level.aggregates.aggregated_path_vals];

finalsale_data_to_use = finalsale_data_to_use(non_empties, :);

split_numerics = arrayfun(@(x) finalsale_data_to_use(x, :), (1:length(finalsale_data_to_use))', 'UniformOutput', false);

finalsale_blocks = cellfun(@(x, y) num2cell(repmat(y, [size(x, 1) 1])), industry_link_indexes, split_numerics, 'UniformOutput', false);
 
linked_consumption_finalsale_data_industry_level = [vertcat(grouped_industry_index_links{:}) vertcat(finalsale_blocks{:}) num2cell(vertcat(scaled_global_industry_threat_intensities{:})) ];

[~, sorted_global_consumption_industry_inds] = sort(cell2mat(linked_consumption_finalsale_data_industry_level(:, end)), 'descend');

linked_consumption_finalsale_data_industry_level = linked_consumption_finalsale_data_industry_level(sorted_global_consumption_industry_inds, :);

linked_consumption_finalsale_data_industry_level = [num2cell((1:length(linked_consumption_finalsale_data_industry_level))') linked_consumption_finalsale_data_industry_level];
dom_consumption_inds = strcmp(linked_consumption_finalsale_data_industry_level(:, 2), linked_consumption_finalsale_data_industry_level(:, 3));

international_linked_consumption_finalsale_data_industry_level = linked_consumption_finalsale_data_industry_level((~dom_consumption_inds), :);

T = cell2table(international_linked_consumption_finalsale_data_industry_level, 'VariableNames', ...
    [{'global_rank', 'consumption_country', 'finalsale_country', 'finalsale_industry'}, analyse_MRIO_params.groups_to_count, {'aggregated_finalsale_threat_intensity', 'finalsale_threat_intensity_proportion', 'attributed_threat_intensity'}]);

writetable(T, '~/GitHub/MRIO_BIO_SATELLITE/international_linked_consumption_finalsale_data_industry_level.txt', 'delimiter', 'tab')

T = cell2table(linked_consumption_finalsale_data_industry_level, 'VariableNames', ...
    [{'global_rank', 'consumption_country', 'finalsale_country', 'finalsale_industry'}, analyse_MRIO_params.groups_to_count, {'aggregated_finalsale_threat_intensity', 'finalsale_threat_intensity_proportion', 'attributed_threat_intensity'}]);

writetable(T, '~/GitHub/MRIO_BIO_SATELLITE/global_linked_consumption_finalsale_data_industry_level.txt', 'delimiter', 'tab')


%%%%%%%%%%%% COUNTRY LEVEL GLOBAL ASSESSMENT


industry_link_indexes = arrayfun(@(x) find(global_aggregated_consumption.industry_level.index_list == x), finalsale_data.industry_level.aggregates.aggregated_paths, 'UniformOutput', false);

non_empties = cellfun('length', industry_link_indexes) > 0;
industry_link_indexes_to_use = industry_link_indexes(non_empties);

finalsale_data_to_use = finalsale_data.industry_level.aggregates.grouped_aggregates(non_empties); %%%%%%%%%%% run analysis on species

%%%%%%%%%%% run analysis on final sale industries
%%%%% finalsale_data_to_use = arrayfun(@(x) {{x}}, finalsale_data.industry_level.aggregates.aggregated_paths(non_empties));

industry_indexes_to_use = finalsale_data.industry_level.aggregates.aggregated_paths(non_empties);

finalsale_blocks = cellfun(@(x, y) repmat({vertcat(y{:})}, [size(x, 1) 1]), industry_link_indexes_to_use, finalsale_data_to_use, 'UniformOutput', false);
finalsale_blocks = vertcat(finalsale_blocks{:});

        
grouped_industry_index_links = cellfun(@(x) [global_aggregated_consumption.industry_level.country_index_list(x) repmat(global_aggregated_consumption.industry_level.index_list(x), [1, 2])], industry_link_indexes_to_use, 'UniformOutput', false);

grouped_industry_index_links = vertcat(grouped_industry_index_links{:});

country_blocks = arrayfun(@(x) finalsale_blocks(grouped_industry_index_links(:, 1) == x), (1:country_num)', 'un', false);
country_blocks = cellfun(@(x) vertcat(x{:}), country_blocks, 'un', false);
country_blocks = cellfun(@(x) unique(x), country_blocks, 'un', false);

country_level_species_counts = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(x), y))), analyse_MRIO_params.groups_to_count), country_blocks, 'un', false);
country_level_species_counts = vertcat(country_level_species_counts{:});

%linked_consumption_finalsale_data_industry_level = [vertcat(grouped_industry_index_links{:}) vertcat(finalsale_blocks{:}) num2cell(vertcat(scaled_global_industry_threat_intensities{:})) ];








% global_aggregated_finalsale_country_list = cellfun(@(x) num2cell(x.industry_level.finalsale), consumption_country_set, 'UniformOutput', false);
% global_aggregated_finalsale_level_threat_intensities = cellfun(@(x) num2cell(x.ranked_threat_proportions_per_industry), consumption_country_set, 'UniformOutput', false);
% 
% industry_link_indexes = arrayfun(@(x) find(global_aggregated_finalsale_country_list == x), finalsale_data.industry_level.aggregated_paths, 'UniformOutput', false);
% 
% non_empties = cellfun('length', industry_link_indexes) > 0;
% industry_link_indexes = industry_link_indexes(non_empties);
% finalsale_data_to_use = cell2mat(finalsale_data.ranked_aggregated_paths(non_empties, 3:5));
% 
% split_numerics = arrayfun(@(x) finalsale_data_to_use(x, :), (1:length(finalsale_data_to_use))', 'UniformOutput', false);
% 
% grouped_industry_index_links = cellfun(@(x) [global_aggregated_consumption.industry_level.country_list(x)...
%                                 IUCN_data_object.industry_characteristics.country_names_list(global_aggregated_finalsale_country_list(x))  ...
%                                IUCN_data_object.industry_characteristics.commodity_classification_list(global_aggregated_finalsale_country_list(x))], industry_link_indexes, 'UniformOutput', false);
%                            
% scaled_global_industry_threat_intensities = cellfun(@(x, y) [global_aggregated_consumption.industry_level.threat_intensities(x) ./ (sum(global_aggregated_consumption.industry_level.threat_intensities(x))) ...
%                                                             global_aggregated_consumption.industry_level.threat_intensities(x) .* (y(3) ./ sum(global_aggregated_consumption.industry_level.threat_intensities(x)))], ...
%                                                           industry_link_indexes, split_numerics, 'UniformOutput', false);
% 
% finalsale_blocks = cellfun(@(x, y) num2cell(repmat(y, [size(x, 1) 1])), industry_link_indexes, split_numerics, 'UniformOutput', false);













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERNATIONAL INDUSTRY ASSESSMENT

unique_country_inds = mat2cell((1:length(unique_countries))', ones(length(unique_countries), 1));
global_consumption.species_level.country_list = cellfun(@(x, y) repmat(y, [size(x.species_level.aggregated_paths, 1) 1]), consumption_country_set, unique_country_inds,'UniformOutput', false);
global_consumption.species_level.country_list = vertcat(global_consumption.species_level.country_list{:});

global_finalsale_production_species_level_index_list = cellfun(@(x) num2cell(x.species_level.aggregated_paths), consumption_country_set, 'UniformOutput', false);
global_finalsale_production_species_level_index_list = cell2mat(vertcat(global_finalsale_production_species_level_index_list{:}));

global_finalsale_production_species_level_threat_intensities = cellfun(@(x) num2cell(x.species_level.aggregated_path_vals), consumption_country_set, 'UniformOutput', false);
global_finalsale_production_species_level_threat_intensities = cell2mat(vertcat(global_finalsale_production_species_level_threat_intensities{:}));

% grouped_industry_index_links = cellfun(@(x) [global_aggregated_consumption_species_level_country_list(x)...
%                                 IUCN_data_object.industry_characteristics.country_index_list(global_aggregated_consumption_species_level_index_list(x))  ...
%                                IUCN_data_object.industry_characteristics.commodity_classification_list(global_aggregated_consumption_species_level_index_list(x))], industry_link_indexes, 'UniformOutput', false);
% 
                          
consumption_finalsale_production_links = [global_consumption.species_level.country_list IUCN_data_object.industry_characteristics.country_index_list(global_finalsale_production_species_level_index_list)];

international_inds = find( sum( abs(diff(consumption_finalsale_production_links(:, [1 3]), 1, 2)), 2) );

net_international_trade_proportion = sum(global_finalsale_production_species_level_threat_intensities(international_inds))/sum(global_finalsale_production_species_level_threat_intensities);

international_consumption_blocks = cellfun(@(x) international_inds( global_consumption.species_level.country_list(international_inds) == x ), unique_country_inds, 'UniformOutput', false);

global_consumption_blocks = cellfun(@(x) find(global_consumption.species_level.country_list == x), unique_country_inds, 'UniformOutput', false);

international_consumption_stats = cellfun(@(x, y) [sum(global_finalsale_production_species_level_threat_intensities(x)) sum(global_finalsale_production_species_level_threat_intensities(y))], ...
    international_consumption_blocks, global_consumption_blocks, 'UniformOutput', false);
international_consumption_stats = vertcat(international_consumption_stats{:});
dom_consumption_trade_stats = international_consumption_stats(:, 2) - international_consumption_stats(:, 1);
global_consumption_trade_stats = [dom_consumption_trade_stats international_consumption_stats international_consumption_stats(:, 1)./international_consumption_stats(:, 2)];
[~, sorted_consumption_inds] = sort(global_consumption_trade_stats(:, 3), 'descend'); 

global_consumption_table = cell2table([unique_countries(sorted_consumption_inds) num2cell(country_level_species_counts(sorted_consumption_inds, :)), num2cell(global_consumption_trade_stats(sorted_consumption_inds, :))], 'VariableNames', ...
    [{'Consumption_Country'}, analyse_MRIO_params.groups_to_count, {'domestic_consumption', 'international_consumption', 'total_consumption', 'international_consumption_proportion'}]);





%%%%%%%%%%%% international production routes

international_production_blocks = cellfun(@(x) international_inds( consumption_finalsale_production_links(international_inds, 3) == x ), unique_country_inds, 'UniformOutput', false);

global_production_blocks = cellfun(@(x) find(consumption_finalsale_production_links(:, 3) == x), unique_country_inds, 'UniformOutput', false);

international_production_stats = cellfun(@(x, y) [sum(global_finalsale_production_species_level_threat_intensities(x)) sum(global_finalsale_production_species_level_threat_intensities(y))], ...
                                                international_production_blocks, global_production_blocks, 'UniformOutput', false);
international_production_stats = vertcat(international_production_stats{:});
dom_production_trade_stats = international_production_stats(:, 2) - international_production_stats(:, 1);
global_production_trade_stats = [dom_production_trade_stats international_production_stats international_production_stats(:, 1)./international_production_stats(:, 2)];
[~, sorted_production_inds] = sort(global_production_trade_stats(:, 2), 'descend');

global_production_table = cell2table([unique_countries(sorted_production_inds) num2cell(global_production_trade_stats(sorted_production_inds, :))], 'VariableNames', ...
    {'Production_Country', 'domestic_production', 'international_production', 'total_production', 'international_production_proportion'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import_export_discriminator = global_consumption_trade_stats(:, 2) - global_production_trade_stats(:, 2);
net_threat_intensity = (global_consumption_trade_stats(:, 3) - global_production_trade_stats(:, 2));

global_net_table = cell2table([unique_countries num2cell([global_consumption_trade_stats(:, [1 2 4]) global_production_trade_stats(:, 2) import_export_discriminator global_consumption_trade_stats(:, 3) net_threat_intensity])], ...
    'VariableNames', {'Consumption_Country', 'domestic_consumption', 'imported_consumption', 'consumption_proportion', 'exported_production', 'net_imported_exported', 'net_consumption','net_threat_intensity'});

global_net_table = global_net_table(sorted_consumption_inds, :);


fid = fopen('~/Github/MRIO_BIO_SATELLITE/developing_countries.txt');
        developing_countries = textscan(fid,'%s', 'HeaderLines', 0, 'delimiter', ';');
fclose(fid);

[~, low_income_inds] = intersect(unique_countries, developing_countries{1});

high_income_inds = (setdiff(1:length(unique_countries), low_income_inds))';

high_income_consumption_low_income_production_indexes = cellfun(@(x) x(ismember(consumption_finalsale_production_links(x, 3), low_income_inds)), international_consumption_blocks(high_income_inds), 'UniformOutput', false);
high_income_consumption_low_income_production_vals = cell2mat(cellfun(@(x) sum(global_finalsale_production_species_level_threat_intensities(x)), high_income_consumption_low_income_production_indexes, 'UniformOutput', false));

low_income_prop = high_income_consumption_low_income_production_vals./global_consumption_trade_stats(high_income_inds, 2);
low_income_stats = [high_income_consumption_low_income_production_vals low_income_prop global_consumption_trade_stats(high_income_inds, 2) ];
[~, sorted_high_incomes] = sort(low_income_stats(:, 3), 'descend');
low_income_table = [unique_countries(high_income_inds) num2cell(low_income_stats)];

T = cell2table(low_income_table(sorted_high_incomes, :), 'VariableNames', {'consumption_country', 'low_income_consumption', 'low_income_proportion', 'net_consumption'});



% high_income_consumption_low_income_sources = cellfun(@(x) ismember(consumption_finalsale_production_links(international_consumption_blocks{x}, 3), low_income_inds), num2cell(high_income_inds'), 'UniformOutput', false);

% international_production_stats = cellfun(@(x, y) sum(global_finalsale_production_species_level_threat_intensities(x(y))), ...
%                                                 international_production_blocks(high_income_inds), high_income_consumption_low_income_sources, 'UniformOutput', false);

                                            
% [~, ranked_global_indexes] = sort(cell2mat(bb(:, 10)), 'descend');
% ranked_global_consumption_industry_paths = bb(ranked_global_indexes, :);

%aa = [IUCN_data_object.industry_characteristics.country_names_list(finalsale_data.industry_level.aggregated_paths) IUCN_data_object.industry_characteristics.commodity_classification_list(finalsale_data.industry_level.aggregated_paths)]
    
%     [aggregated_paths, ~, path_indexes] = unique(paths_to_aggregate, 'rows', 'stable');
%     
%     aggregated_path_vals = accumarray(path_indexes, vals_to_aggregate);
%     
%     grouped_path_indexes = accumarray(path_indexes, find(path_indexes), [], @(rows){rows});  
%     [aggregated_path_vals, sorted_indexes] = sort(aggregated_path_vals, 'descend');
%     aggregated_paths = aggregated_paths(sorted_indexes, :);
%     grouped_path_indexes = grouped_path_indexes(sorted_indexes);
% 
%     sorted_grouped_path_indexes = cellfun(@(x) x(sort_indexes(vals_to_aggregate(x))), grouped_path_indexes, 'un', false);
%     grouped_path_vals = cellfun(@(x) vals_to_aggregate(x), sorted_grouped_path_indexes, 'un', false);
%     grouped_aggregates = cellfun(@(x) objects_to_aggregate(x, :), sorted_grouped_path_indexes, 'un', false);
    





%writetable(name_data, '~/GitHub/MRIO_BIO_SATELLITE/name_data.txt', 'delimiter', 'tab')
%writetable(t, '~/GitHub/MRIO_BIO_SATELLITE/global_rank_list_with_consumption_country.txt', 'delimiter', 'tab')
 