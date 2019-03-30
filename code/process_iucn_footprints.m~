function consumption_footprint_characteristics = process_iucn_footprints(country_of_interest, impact_assessment_level, industry_inputs, footprint_input_objects, analyse_footprint_params, build_footprint_params)
    
    [unique_countries, country_indexes_to_use] = setdiff(industry_inputs.industry_characteristics.unique_countries, analyse_footprint_params.countries_to_exclude, 'stable');    
    
%    disp('processing footprints at final sale level...')   
%     processed_finalsale_data = cellfun(@(x) process_footprint_at_finalsale_level(industry_characteristics, analyse_footprint_params, footprint_filename_prefix, species_characteristics, ...
%                                                                                  footprint_input_objects.satellite_collapse_groups{x}, analyse_footprint_params.species_group_weights(x)), ...
%                                                                                  num2cell(1:numel(footprint_objects.species_group_names))', 'un', false);
                           
    consumption_footprint_characteristics = cellfun(@(x, y) build_consumption_level_footprints(country_of_interest, impact_assessment_level, analyse_footprint_params, build_footprint_params, industry_inputs, unique_countries, ...
                                                                            country_indexes_to_use, footprint_input_objects, x), footprint_input_objects.species_group_names, 'un', false);
             
    consumption_footprint_characteristics = cell2struct(consumption_footprint_characteristics, footprint_input_objects.species_group_names, 1);
    
    consumption_footprint_characteristics.aggregated.country_scale.consumption_characteristics = aggregate_grouped_data(consumption_footprint_characteristics, unique_countries, footprint_input_objects, ...
                                                                                                        'country_scale', 'consumption_characteristics');
    consumption_footprint_characteristics.aggregated.country_scale.net_trade_characteristics = aggregate_grouped_data(consumption_footprint_characteristics, unique_countries, footprint_input_objects, ...
                                                                                                       'country_scale', 'net_trade_characteristics');
    consumption_footprint_characteristics.aggregated.country_scale.production_characteristics = aggregate_grouped_data(consumption_footprint_characteristics, unique_countries, footprint_input_objects, ...
                                                                                                        'country_scale', 'production_characteristics');

end


function aggregated_data = aggregate_grouped_data(trade_characteristics, unique_countries, footprint_input_objects, trade_scale, field_to_aggregate)
    
    [blocknames, data_to_aggregate] = cellfun(@(x, y) balance_data_block(trade_characteristics.(x).(trade_scale).(field_to_aggregate), ...
                                                        trade_characteristics.(x).included_countries, unique_countries, numel(y)), ... 
                                                        footprint_input_objects.species_group_names, footprint_input_objects.satellite_collapse_groups, 'un', false);
                                  
    aggregated_data = sum(cat(3, data_to_aggregate{:}),3);
    aggregated_data = cell2table([unique_countries num2cell(aggregated_data)], 'VariableNames', [{'Country'} blocknames{1}]);
    
end


function [blocknames, sorted_data] = balance_data_block(current_group, current_country_set, unique_countries, species_num)
    
    blocknames = current_group.data_blocknames;
    sorted_data = zeros(numel(unique_countries), numel(blocknames));
    [~, template_indexes, sorted_data_indexes] = intersect(unique_countries, current_country_set, 'stable');
    
    sorted_data(template_indexes, :) = current_group.data_block(sorted_data_indexes, :)/species_num;

%     norm_const = sum(sorted_data(:, strcmp(current_group.data_blocknames, 'total')));
%     
%     groups_to_balance = ~ismember(current_group.data_blocknames, species_group_names);
%     
%     sorted_data(:, groups_to_balance) = current_weight * sorted_data(:, groups_to_balance)/norm_const;
    
end


function consumption_footprint_data = build_consumption_level_footprints(country_of_interest, impact_assessment_level, analyse_footprint_params, build_footprint_params, industry_inputs, ...
                                                                                unique_countries, country_indexes_to_use, footprint_input_objects, species_group_to_use)
         
    display(['calculating footprints for ', species_group_to_use ' species'])
    consumption_footprint_data = struct();
    
    grouped_finalsale_data = load_consumption_footprints(country_of_interest, ...
                                                         'consumption', ...
                                                          analyse_footprint_params, ...
                                                          build_footprint_params, ...
                                                          impact_assessment_level, ...
                                                          industry_inputs.industry_characteristics, ...
                                                          country_indexes_to_use, ...
                                                          footprint_input_objects, ...
                                                          species_group_to_use);
    
    species_present = find(cellfun('length', grouped_finalsale_data));  
    
    if ~isempty(species_present)
        
        country_indexes_to_use = country_indexes_to_use(species_present);

        consumption_footprint_data.included_countries = unique_countries(species_present);
        
        grouped_finalsale_data = grouped_finalsale_data(species_present);   
        
        consumption_footprint_data.sector_to_sector_scale = expand_sector_to_sector_data(grouped_finalsale_data, ...
                                                                                                  country_indexes_to_use, ...
                                                                                                  industry_inputs.industry_characteristics, ...
                                                                                                  impact_assessment_level, ...
                                                                                                  country_of_interest);
                                                                                              
        consumption_footprint_data.aggregated_sector_scale = expand_aggregated_sector_scale_data(grouped_finalsale_data, ...
                                                                                                 country_indexes_to_use, ... 
                                                                                                 industry_inputs.industry_characteristics, ...
                                                                                                 impact_assessment_level,... 
                                                                                                 country_of_interest);
                                                                                             
        consumption_footprint_data.country_scale.consumption_characteristics = assess_international_trade_characteristics('consumption', ...
                                                                                                                          consumption_footprint_data,...
                                                                                                                          consumption_footprint_data.included_countries, ...
                                                                                                                          country_indexes_to_use); 
                                                                                                                      
        consumption_footprint_data.country_scale.production_characteristics = assess_international_trade_characteristics('production', ...
                                                                                                               consumption_footprint_data, ...
                                                                                                               consumption_footprint_data.included_countries, ...
                                                                                                               country_indexes_to_use);
                                                                                       
        consumption_footprint_data.country_scale.import_characteristics = build_import_trade_characteristics(consumption_footprint_data, ...
                                                                                                             industry_inputs.low_income_country_names, ...   
                                                                                                             consumption_footprint_data.included_countries, ....
                                                                                                             country_indexes_to_use);
    
        consumption_footprint_data.country_scale.net_trade_characteristics = build_country_scale_net_trade_table(consumption_footprint_data, ...
                                                                                                                 consumption_footprint_data.included_countries);     
    
        consumption_footprint_data.global_scale = build_global_scale_characteristics(consumption_footprint_data);
        
    end
    
end

        
function aggregated_outputs = aggregate_footprint_at_finalsale_routines(current_footprint, analyse_footprint_params, industry_characteristics)
                     
        aggregated_outputs.sector_to_sector_scale = aggregate_and_sort_paths(analyse_footprint_params.sort_data, ...
                                                      current_footprint.vals, ...
                                                      [current_footprint.finalsale current_footprint.production], ...
                                                      current_footprint.species);
   

    if ~isempty(fieldnames(aggregated_outputs.sector_to_sector_scale))
        
        aggregated_outputs.aggregated_sector_scale = aggregate_and_sort_paths(analyse_footprint_params.sort_data, ...
                                                                              aggregated_outputs.sector_to_sector_scale.aggregated_path_vals, ...
                                                                              aggregated_outputs.sector_to_sector_scale.aggregated_paths(:, 1), ...
                                                                              aggregated_outputs.sector_to_sector_scale.grouped_aggregates);
   
        aggregated_outputs.aggregated_sector_scale.grouped_sectors = cellfun(@(x) aggregated_outputs.sector_to_sector_scale.aggregated_paths(x, 2), ...
                                                                            aggregated_outputs.aggregated_sector_scale.grouped_path_indexes, 'un', false);
        
        aggregated_outputs.country_scale.production = aggregate_and_sort_paths(analyse_footprint_params.sort_data, ...
                                                                    aggregated_outputs.sector_to_sector_scale.aggregated_path_vals, ...
                                                                    industry_characteristics.country_index_list(aggregated_outputs.sector_to_sector_scale.aggregated_paths(:, 2)), ...
                                                                    aggregated_outputs.sector_to_sector_scale.grouped_aggregates);
  
        aggregated_outputs.country_scale.finalsale = aggregate_and_sort_paths(analyse_footprint_params.sort_data, ...
                                                                   aggregated_outputs.sector_to_sector_scale.aggregated_path_vals, ...
                                                                   industry_characteristics.country_index_list(aggregated_outputs.sector_to_sector_scale.aggregated_paths(:, 1)), ...
                                                                   aggregated_outputs.sector_to_sector_scale.grouped_aggregates);   
    end
    
end







function species_group_counts = count_species_groups(grouped_species_aggregates, species_characteristics, groups_to_count)

    species_group_counts = cell(numel(groups_to_count),1);
    species_group_counts = cell2struct(species_group_counts, groups_to_count);
    for current_counter = 1:length(groups_to_count)
        species_to_use = cellfun(@(x) unique(vertcat(x{:})), grouped_species_aggregates, 'un', false);
        species_group_counts.(groups_to_count{current_counter}) = cell2mat(cellfun(@(x) length(find(strcmp(species_characteristics.species_kingdom(x), groups_to_count{current_counter}))), ...
                                                         species_to_use, 'un', false));
    end
    
end


function [outputs] = aggregate_and_sort_paths(sort_data, vals_to_aggregate, paths_to_aggregate, objects_to_aggregate)
    
    outputs = struct();
    [aggregated_paths, ~, path_indexes] = unique(paths_to_aggregate, 'rows', 'stable');
    
    if ~isempty(path_indexes)
        
        outputs.aggregated_paths = aggregated_paths;
        outputs.aggregated_path_vals = accumarray(path_indexes, vals_to_aggregate);
    
        outputs.grouped_path_indexes = accumarray(path_indexes, find(path_indexes), [], @(rows){rows});
        
        if sort_data == true
            [outputs.aggregated_path_vals, sorted_indexes] = sort(outputs.aggregated_path_vals, 'descend');
            outputs.aggregated_paths = outputs.aggregated_paths(sorted_indexes, :);
            outputs.grouped_path_indexes = outputs.grouped_path_indexes(sorted_indexes);
            outputs.grouped_path_indexes = cellfun(@(x) x(sort_indexes(vals_to_aggregate(x))), outputs.grouped_path_indexes, 'un', false);
        end
    
        outputs.grouped_path_vals = cellfun(@(x) vals_to_aggregate(x), outputs.grouped_path_indexes, 'un', false);
        outputs.grouped_aggregates = cellfun(@(x) objects_to_aggregate(x, :), outputs.grouped_path_indexes, 'un', false);

    end
   
   
end

function [sorted_inds] = sort_indexes(arr)
    [~, sorted_inds] = sort(arr, 'descend');
end



function current_footprint = select_footprint_subset(current_footprint, country_of_interest, impact_assessment_level, species_group_to_use, industry_characteristics, threshold_level)
    
    if (strcmp(species_group_to_use, 'all'))
        species_group_to_use = 1:numel(current_footprint.species);
    end
        
%     current_footprint = vertcat(current_footprint{[species_group_to_use]});
%     current_footprint = structfun(@(x) vertcat(x{:}), current_footprint);
%     
%     current_footprint = vertcat(cf{[1 2]});
%     current_footprint = structfun(@(x) vertcat(x), current_footprint);
    
    current_footprint = structfun(@(x) vertcat(x{species_group_to_use}), current_footprint, 'un', false);
  
    if ~isempty(current_footprint.species)
       species_group_to_use = current_footprint.vals > threshold_level;
       current_footprint = structfun(@(x) x(species_group_to_use), current_footprint, 'un', false);
    end
    
    if ( ~isempty(current_footprint.species) && ~strcmp(country_of_interest, 'global') && ~strcmp(impact_assessment_level, 'consumption'))
        
        industries_to_use = find(ismember(industry_characteristics.country_names_list, country_of_interest));
        
        if strcmp(impact_assessment_level, 'production')
            current_industry_set = ismember(current_footprint.production, industries_to_use);
        elseif strcmp(impact_assessment_level, 'finalsale')
            current_industry_set = ismember(current_footprint.finalsale, industries_to_use);
        end
        
        current_footprint = structfun(@(x) x(current_industry_set), current_footprint, 'un', false);
    
    end
    
end


function [species_threat_num] = build_species_threat_num(iucn_species_names, unique_iucn_species)
    species_threat_num = zeros(numel(unique_iucn_species), 1);
    
    for i = 1:numel(unique_iucn_species)
        species_threat_num(i) = sum(strcmp(iucn_species_names, unique_iucn_species(i)));
    end

end


function global_scale_characteristics = build_global_scale_characteristics(trade_characteristics)
    display('....at global scale ')
    
    global_scale_characteristics = struct();
    global_scale_characteristics.trade_stats.total_global_trade = sum(trade_characteristics.sector_to_sector_scale.finalsale_to_production_threat_intensities);
    global_scale_characteristics.trade_stats.total_international_trade = sum(trade_characteristics.sector_to_sector_scale.finalsale_to_production_threat_intensities(trade_characteristics.sector_to_sector_scale.international_indexes));
    global_scale_characteristics.trade_stats.proportion_international_trade = global_scale_characteristics.trade_stats.total_international_trade/global_scale_characteristics.trade_stats.total_global_trade;
    
    global_scale_characteristics.international_trade.low_income = calc_global_stats(trade_characteristics.country_scale.import_characteristics.low_income, ...
                                                                                    global_scale_characteristics.trade_stats.total_international_trade);
     global_scale_characteristics.international_trade.high_income = calc_global_stats(trade_characteristics.country_scale.import_characteristics.high_income, ...
                                                                                    global_scale_characteristics.trade_stats.total_international_trade);
    
end

function global_scale_characteristics = calc_global_stats(trade_object, total_international_trade)
    
    global_scale_characteristics = struct();
    global_scale_characteristics.total = sum(trade_object.imported_total);
    global_scale_characteristics.total_proportion_of_international = global_scale_characteristics.total./(total_international_trade + 1e-10);
    global_scale_characteristics.proportion_imported_by_high = sum(trade_object.imported_by_high)./(global_scale_characteristics.total + 1e-10);
    global_scale_characteristics.proportion_imported_by_low = sum(trade_object.imported_by_low)./(global_scale_characteristics.total + 1e-10);
    
end


function country_scale_consumption = load_consumption_footprints(country_of_interest, footprint_level, analyse_footprint_params, build_footprint_params, impact_assessment_level, industry_characteristics, country_indexes_to_use, ...
                                                                 footprint_input_objects, species_group_to_use)
                                                                        
    country_num = numel(country_indexes_to_use);
    country_scale_consumption = cell(country_num, 1);
    [~, species_group_numerical] = intersect(footprint_input_objects.species_group_names, species_group_to_use);
    
    for country_index = 1:country_num
        
        current_footprint_filename =  [build_footprint_params.footprint_filename_prefix footprint_level '_' industry_characteristics.unique_country_codes{country_indexes_to_use(country_index)}  ...
                                  build_footprint_params.footprint_filename_suffix];
        
        current_footprint = load(current_footprint_filename);
                              
        obj_names = fieldnames(current_footprint);
        current_footprint = current_footprint.(obj_names{1});
        
        current_footprint = select_footprint_subset(current_footprint, ...
                                                    country_of_interest, ...
                                                    impact_assessment_level, ...
                                                    species_group_numerical, ...
                                                    industry_characteristics, ...
                                                    analyse_footprint_params.local_thresh_level); 
                                                         
        if (length(current_footprint.species) > 1)
        	country_scale_consumption{country_index} = aggregate_footprint_at_finalsale_routines(current_footprint, analyse_footprint_params, industry_characteristics);
        end

    end
    
end



function aggregated_sector_scale = expand_aggregated_sector_scale_data(consumption_country_set, country_indexes_to_use, industry_characteristics, impact_assessment_level, country_of_interest)
    
    display('....at aggregated sector scale ')
    aggregated_sector_scale = struct();
    aggregated_sector_scale.consumption_country_index_list = cellfun(@(x, y) repmat(y, [size(x.aggregated_sector_scale.aggregated_paths, 1) 1]), consumption_country_set, num2cell(country_indexes_to_use), 'un', false);
    aggregated_sector_scale.finalsale_sector_list = cellfun(@(x) x.aggregated_sector_scale.aggregated_paths, consumption_country_set, 'un', false);
    aggregated_sector_scale.threat_intensities = cellfun(@(x) x.aggregated_sector_scale.aggregated_path_vals, consumption_country_set, 'un', false);
    
    aggregated_sector_scale = structfun(@(x) vertcat(x{:}), aggregated_sector_scale, 'un', false);
    
    if strcmp(impact_assessment_level, 'consumption') && ~strcmp(country_of_interest, 'global')
        consumption_countries_use = find(ismember(industry_characteristics.unique_countries, country_of_interest));
        current_industry_set = ismember(aggregated_sector_scale.consumption_country_index_list, consumption_countries_use);
        aggregated_sector_scale = structfun(@(x) x(current_industry_set), aggregated_sector_scale, 'un', false);
    end
    
    data_block = [industry_characteristics.unique_countries(aggregated_sector_scale.consumption_country_index_list), ...
                                     industry_characteristics.country_names_list(aggregated_sector_scale.finalsale_sector_list), ...
                                     industry_characteristics.commodity_classification_list(aggregated_sector_scale.finalsale_sector_list), ...
                                     num2cell(aggregated_sector_scale.threat_intensities)];
                                 
    aggregated_sector_scale.table = cell2table(data_block, 'VariableNames', {'consumption_country', 'finalsale_country', 'finalsale_industry', 'threat_intensity'});
    
end


function sector_to_sector_scale = expand_sector_to_sector_data(consumption_country_set, country_indexes_to_use, industry_characteristics, impact_assessment_level, country_of_interest)
    
    display('....at sector scale ')
    sector_to_sector_scale = struct();
    
    sector_to_sector_scale.finalsale_to_production_sectors = cellfun(@(x) x.sector_to_sector_scale.aggregated_paths, consumption_country_set, 'un', false); 
    sector_to_sector_scale.finalsale_to_production_threat_intensities = cellfun(@(x) x.sector_to_sector_scale.aggregated_path_vals, consumption_country_set, 'un', false);
    sector_to_sector_scale.consumption_country_index_list = cellfun(@(x, y) repmat(y, [size(x, 1) 1]), sector_to_sector_scale.finalsale_to_production_sectors, num2cell(country_indexes_to_use), 'un', false);

    sector_to_sector_scale = structfun(@(x) vertcat(x{:}), sector_to_sector_scale, 'un', false);
        
    sector_to_sector_scale.finalsale_country_list = industry_characteristics.country_index_list(sector_to_sector_scale.finalsale_to_production_sectors(:, 1));
    sector_to_sector_scale.production_country_list = industry_characteristics.country_index_list(sector_to_sector_scale.finalsale_to_production_sectors(:, 2));
    
     if strcmp(impact_assessment_level, 'consumption') && ~strcmp(country_of_interest, 'global')
        consumption_countries_to_use = find(ismember(industry_characteristics.unique_countries, country_of_interest));
        current_industry_set = ismember(sector_to_sector_scale.consumption_country_index_list, consumption_countries_to_use);
        sector_to_sector_scale = structfun(@(x) x(current_industry_set), sector_to_sector_scale, 'un', false);
     end
    
    sector_to_sector_scale.international_indexes = sum( abs(diff([sector_to_sector_scale.consumption_country_index_list ...
                                                                  sector_to_sector_scale.finalsale_country_list ...
                                                                  sector_to_sector_scale.production_country_list], 1, 2)), 2) > 0;
     
    data_block = [industry_characteristics.unique_countries(sector_to_sector_scale.consumption_country_index_list), ...
                                     industry_characteristics.country_names_list(sector_to_sector_scale.finalsale_sector_list), ...
                                     industry_characteristics.commodity_classification_list(sector_to_sector_scale.finalsale_sector_list), ...
                                     num2cell(sector_to_sector_scale.threat_intensities)];
                                 
    sector_to_sector_scale.table = cell2table(data_block, 'VariableNames', {'consumption_country', 'finalsale_country', 'finalsale_industry', 'threat_intensity'});
                                                             
end


function country_scale_net_trade_characteristics = build_country_scale_net_trade_table(trade_characteristics, included_countries)
 
    country_scale_net_trade_characteristics = struct();
    country_scale_net_trade_characteristics.import_export_discriminator = trade_characteristics.country_scale.consumption_characteristics.partitioned_threat_intensities.international ...
                                                                            - trade_characteristics.country_scale.production_characteristics.partitioned_threat_intensities.international;   
                                            
    country_scale_net_trade_characteristics.net_threat_intensity = trade_characteristics.country_scale.consumption_characteristics.partitioned_threat_intensities.global ...
                                            - trade_characteristics.country_scale.production_characteristics.partitioned_threat_intensities.international;
                                                 
    country_scale_net_trade_characteristics.data_block = [trade_characteristics.country_scale.consumption_characteristics.data_block(:, 1:(end - 1))...
                                                     trade_characteristics.country_scale.import_characteristics.low_income.imported_total...
                                                     trade_characteristics.country_scale.production_characteristics.partitioned_threat_intensities.international...
                                                     trade_characteristics.country_scale.consumption_characteristics.data_block(:, end)...
                                                     country_scale_net_trade_characteristics.net_threat_intensity];
                                                 
    country_scale_net_trade_characteristics.data_blocknames = [trade_characteristics.country_scale.consumption_characteristics.data_blocknames(1:(end - 1)) ...
                                                          {'low_income_imports', 'exports', 'gross'} trade_characteristics.country_scale.consumption_characteristics.data_blocknames(end)];

    country_scale_net_trade_characteristics.table = cell2table( [included_countries num2cell([country_scale_net_trade_characteristics.data_block ...
                                                     trade_characteristics.country_scale.import_characteristics.low_income.imported_proportion...
                                                     trade_characteristics.country_scale.production_characteristics.partitioned_threat_intensities.international_trade_proportion...
                                                     trade_characteristics.country_scale.consumption_characteristics.partitioned_threat_intensities.international_trade_proportion])],...
                                                     'VariableNames', [{'Country'} country_scale_net_trade_characteristics.data_blocknames {'low_income_imported_proportion', 'proportion_exported', 'international_consumption_proportion'}] );
    
end


function country_scale_trade_characteristics = assess_international_trade_characteristics(assess_type, consumption_footprint_data, included_countries, country_indexes_to_use)
    display('....at country scale ')
    if strcmp(assess_type, 'consumption')
        country_index_list = consumption_footprint_data.sector_to_sector_scale.consumption_country_index_list;   
    else
       country_index_list = consumption_footprint_data.sector_to_sector_scale.production_country_list;
    end  

    country_blocks = cellfun(@(x) country_index_list == x, num2cell(country_indexes_to_use), 'un', false);
    
    country_scale_trade_characteristics.partitioned_indexes.global = cellfun(@(x) find(x), country_blocks, 'un', false);
    country_scale_trade_characteristics.partitioned_indexes.international = cellfun(@(x) find(x & consumption_footprint_data.sector_to_sector_scale.international_indexes), country_blocks, 'un', false);
    country_scale_trade_characteristics.partitioned_indexes.domestic = cellfun(@(x) find(x & ~consumption_footprint_data.sector_to_sector_scale.international_indexes), country_blocks, 'un', false);
    
    country_scale_trade_characteristics.partitioned_threat_intensities = structfun(@(x) cellfun(@(y) sum(consumption_footprint_data.sector_to_sector_scale.finalsale_to_production_threat_intensities(y)), x), ...
                                                                                        country_scale_trade_characteristics.partitioned_indexes, 'un', false);
    
    
    country_scale_trade_characteristics.partitioned_threat_intensities.international_trade_proportion = country_scale_trade_characteristics.partitioned_threat_intensities.international ./ ...
                                                                                                        (country_scale_trade_characteristics.partitioned_threat_intensities.global + 1e-10);
        
    country_scale_trade_characteristics.data_block = [country_scale_trade_characteristics.partitioned_threat_intensities.domestic ...
                                                      country_scale_trade_characteristics.partitioned_threat_intensities.international ...
                                                      country_scale_trade_characteristics.partitioned_threat_intensities.global];
    
    country_scale_trade_characteristics.data_blocknames = {'domestic', 'international', 'net'};
    country_scale_trade_characteristics.table = cell2table([included_countries num2cell(country_scale_trade_characteristics.data_block) ...
                                                            num2cell(country_scale_trade_characteristics.partitioned_threat_intensities.international_trade_proportion)],...
                                                            'VariableNames', [{'Country'}, country_scale_trade_characteristics.data_blocknames  {'international_proportion'}]);
end


function import_characteristics = build_import_trade_characteristics(trade_characteristics, low_income_country_names, included_countries, country_indexes_to_use)
    
    import_characteristics = struct();
    
    [~, low_income_internal_indexes] = intersect(included_countries, low_income_country_names{1});
    
    high_income_internal_indexes = setdiff(1:length(included_countries), low_income_internal_indexes);
    
    import_characteristics.low_income = assess_imports_by_economy_level('low', ...
                                                                        trade_characteristics, ...
                                                                        trade_characteristics.country_scale.consumption_characteristics.partitioned_indexes.international, ...
                                                                        low_income_internal_indexes, ...
                                                                        high_income_internal_indexes, ...
                                                                        country_indexes_to_use, ....
                                                                        included_countries);
                                                                    
    import_characteristics.high_income = assess_imports_by_economy_level('high', ...
                                                                         trade_characteristics, ...
                                                                         trade_characteristics.country_scale.consumption_characteristics.partitioned_indexes.international, ... 
                                                                         low_income_internal_indexes, ...
                                                                         high_income_internal_indexes, ...
                                                                         country_indexes_to_use, ...
                                                                         included_countries);
    
end


function import_characteristics = assess_imports_by_economy_level(income_type, trade_characteristics, international_trade_indexes, low_income_internal_indexes, high_income_internal_indexes, country_indexes_to_use, unique_countries)
    
    import_characteristics = struct();
    
    if strcmp(income_type, 'low')
        internal_indexes_to_use = low_income_internal_indexes;
    else internal_indexes_to_use = high_income_internal_indexes;
    end
    
    import_characteristics.country_names = unique_countries(internal_indexes_to_use);
    import_characteristics.import_country_indexes = country_indexes_to_use(internal_indexes_to_use);
    
    import_indexes_tmp = cellfun(@(x) ismember(trade_characteristics.sector_to_sector_scale.production_country_list(x), import_characteristics.import_country_indexes), ...
                                                               international_trade_indexes, 'un', false);
                                                           
    import_characteristics.import_indexes = cellfun(@(x, y) x(y), international_trade_indexes, import_indexes_tmp, 'un', false);
    
    import_characteristics.imported_total = cellfun(@(x) sum(trade_characteristics.sector_to_sector_scale.finalsale_to_production_threat_intensities(x)), import_characteristics.import_indexes);
     
    import_characteristics.imported_by_low = import_characteristics.imported_total(low_income_internal_indexes);
    import_characteristics.imported_by_high = import_characteristics.imported_total(high_income_internal_indexes);

    import_characteristics.imported_proportion = import_characteristics.imported_total ./ ...
                                                 (trade_characteristics.country_scale.consumption_characteristics.partitioned_threat_intensities.international + 1e-10);
    
    import_characteristics.data_block = [import_characteristics.imported_total trade_characteristics.country_scale.consumption_characteristics.partitioned_threat_intensities.global];
    import_characteristics.data_blocknames = {'imports', 'net_consumption'}; 
    import_characteristics.table = cell2table([unique_countries num2cell(import_characteristics.data_block) num2cell(import_characteristics.imported_proportion)], ...
                                               'VariableNames', [{'Country'}, import_characteristics.data_blocknames {'proportional_imports'}]);
end


function species_counts = build_country_scale_species_counts(merged_finalsale_data, finalsale_sector_list, consumption_country_index_list, finalsale_data, country_indexes_to_use, species_groups)

    industry_link_indexes = arrayfun(@(x) find(finalsale_sector_list == x), finalsale_data.aggregated_sector_scale.aggregated_paths, 'un', false);
    
    non_empties = cellfun('length', industry_link_indexes) > 0;
    industry_link_indexes_to_use = industry_link_indexes(non_empties);

    %%%%%%%%%%% run species count analysis
    finalsale_data_to_use = finalsale_data.aggregated_sector_scale.grouped_aggregates(non_empties); 

    finalsale_blocks = cellfun(@(x, y) repmat({vertcat(y{:})}, [size(x, 1) 1]), industry_link_indexes_to_use, finalsale_data_to_use, 'un', false);
    finalsale_blocks = vertcat(finalsale_blocks{:});
     
    grouped_industry_index_links = cellfun(@(x) [consumption_country_index_list(x) repmat(finalsale_sector_list(x), [1, 2])], industry_link_indexes_to_use, 'un', false);

    grouped_industry_index_links = vertcat(grouped_industry_index_links{:});

    country_species_blocks = arrayfun(@(x) finalsale_blocks(grouped_industry_index_links(:, 1) == x), country_indexes_to_use, 'un', false);
    country_species_blocks = cellfun(@(x) vertcat(x{:}), country_species_blocks, 'un', false);
    country_species_blocks = cellfun(@(x) unique(x), country_species_blocks, 'un', false);

    species_counts = cellfun(@(x) (cellfun(@(y) length(find(ismember(x, y))), species_groups))', country_species_blocks, 'un', false);
    species_counts = vertcat(species_counts{:});
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ROUTINES TO MATCH FINALSALE TO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CONSUMPTION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%         consumption_footprint_data.country_scale.species_counts = build_country_scale_species_counts(consumption_footprint_data.aggregated_sector_scale, ...
%                                                                                       consumption_footprint_data.aggregated_sector_scale.finalsale_sector_list, ... 
%                                                                                       consumption_footprint_data.aggregated_sector_scale.consumption_country_index_list, ... 
%                                                                                       processed_finalsale_data, ...
%                                                                                       country_indexes_to_use, ...
%                                                                                       species_groups);
%
%         consumption_finalsale_matches = match_consumption_to_finalsale(consumption_footprint_data, processed_finalsale_data, country_indexes_to_use);
% 
%         merged_data = struct();
%         merged_data.included_countries = unique_countries(species_present);
% 
%         merged_data.aggregated_sector_scale = merge_consumption_data_with_finalsale(consumption_finalsale_matches, ...
%                                                                                     processed_finalsale_data, ...   
%                                                                                     analyse_footprint_params, ...
%                                                                                     industry_inputs.industry_characteristics);
%
%             merged_data.country_scale.species_counts = build_country_scale_species_counts(consumption_footprint_data.aggregated_sector_scale.finalsale_sector_list, ... 
%                                                                                       consumption_footprint_data.aggregated_sector_scale.consumption_country_index_list, ... 
%                                                                                       processed_finalsale_data, ...
%                                                                                       country_indexes_to_use, ...
%                                                                                       species_groups);
%
%         merged_data.country_scale.consumption_characteristics = assess_international_trade_characteristics('consumption', ...
%                                                                                                             consumption_footprint_data,...
%                                                                                                                 merged_data.included_countries, ...
%                                                                                                                 country_indexes_to_use, ...
%                                                                                                                 merged_data.country_scale.species_counts, ...
%                                                                                                                 species_group_names);
%     
% 
% 
% 
%

% 
% function merged_finalsale_data = build_consumption_finalsale_production_blocks(consumption_country_grouped_by_commodity, matched_finalsale_data, matched_production_data)
% 
%     finalsale_matched_to_production_commodities = cellfun(@(x, y) repmat(x, [numel(y), 1]), num2cell(matched_finalsale_data), matched_production_data, 'un', false);
%     
%     merged_finalsale_data.finalsale = cellfun(@(x, y) repmat({x}, [numel(y), 1]), finalsale_matched_to_production_commodities, consumption_country_grouped_by_commodity, 'un', false);
%     merged_finalsale_data.production = cellfun(@(x, y) repmat({x}, [numel(y), 1]), matched_production_data, consumption_country_grouped_by_commodity, 'un', false);
%     merged_finalsale_data.consumption_country = cellfun(@(x, y) arrayfun(@(z) repmat(z, [size(y, 1) 1]), x, 'un', false), ...
%                                                                consumption_country_grouped_by_commodity, matched_production_data, 'un', false);
%            
%     merged_finalsale_data = structfun(@(x) vertcat(x{:}), merged_finalsale_data, 'un', false);
%     merged_finalsale_data = structfun(@(x) vertcat(x{:}), merged_finalsale_data, 'un', false);
%       
% end
% 
% 
% function merged_finalsale_data = build_scaled_threat_intensities(consumption_country_grouped_by_commodity, finalsale_data, ...
%                                                                  matched_consumption_finalsale_commodity_indexes, matched_finalsale_threat_intensities)
%     
%      [~, ~, matched_consumption_finalsale_commodity_indexes] = intersect(consumption_country_grouped_by_commodity.aggregated_paths, ...
%                                                                             finalsale_data.aggregated_sector_scale.aggregated_paths, 'stable');
%     
%     sorted_finalsale_sector_data = [industry_characteristics.country_names_list(finalsale_data.aggregated_sector_scale.aggregated_paths(sorted_indexes)),...
%                                     industry_characteristics.commodity_classification_list(finalsale_data.aggregated_sector_scale.aggregated_paths(sorted_indexes)) ...
%                                     num2cell(finalsale_data.aggregated_sector_scale.aggregated_path_vals(sorted_indexes))]  ;    
% 
%     [~, sorted_indexes] = sort(consumption_country_grouped_by_commodity.aggregated_path_vals, 'descend');
%     sorted_consumption_data = [industry_characteristics.country_names_list(consumption_country_grouped_by_commodity.aggregated_paths(sorted_indexes)),...
%                                     industry_characteristics.commodity_classification_list(consumption_country_grouped_by_commodity.aggregated_paths(sorted_indexes)), ...
%                                     num2cell(consumption_country_grouped_by_commodity.aggregated_path_vals(sorted_indexes))]  ;    
%     %consumption_country_grouped_by_commodity_intensities = consumption_country_grouped_by_commodity.grouped_path_vals;
%     
%     %scaled_threat_proportion = cellfun(@(x) x./(sum(x)), consumption_country_grouped_by_commodity_intensities, 'un', false);  
% 
%     %scaled_threat_intensity = cellfun(@(x, y) arrayfun(@(z) z.*y, x, 'un', false), scaled_threat_proportion, matched_finalsale_threat_intensities, 'un', false);
%    
%     expanded_finalsale_paths = cellfun(@(x, y) repmat(x, [size(y), 1]), num2cell(consumption_country_grouped_by_commodity.aggregated_paths), ...
%                                         consumption_country_grouped_by_commodity.grouped_aggregates, 'un', false);
%                                     
%     expanded_finalsale_vals = cellfun(@(x, y) repmat(x, [size(y), 1]), num2cell(consumption_country_grouped_by_commodity.aggregated_path_vals), ...
%                                         consumption_country_grouped_by_commodity.grouped_aggregates, 'un', false);  
%     
%     test2 = [vertcat(expanded_finalsale_vals{:}) vertcat(expanded_finalsale_paths{:})] ;
%     
%     cellfun(@(x, y) repmat(x, [size(y), 1]), num2cell(finalsale_data.aggregated_sector_scale.aggregated_paths), ...
%                                         consumption_country_grouped_by_commodity.grouped_aggregates, 'un', false);
%                                     
%     %test_data = [vertcat(consumption_country_grouped_by_commodity.grouped_aggregates{:}) vertcat(consumption_country_grouped_by_commodity.grouped_path_vals{:}) ]
%     
%     %us_inds = find(test_data(:, 1) == 181);
%     
%     consumption_country_grouped_by_commodity_intensities = consumption_country_grouped_by_commodity.grouped_path_vals;
%     merged_finalsale_data.scaled_threat_proportion = cellfun(@(x) x./(sum(x) + 1e-10), consumption_country_grouped_by_commodity_intensities, 'un', false);  
% 
%     scaled_threat_intensity = cellfun(@(x, y) arrayfun(@(z) z.*y, x, 'un', false), merged_finalsale_data.scaled_threat_proportion, num2cell(matched_finalsale_threat_intensities), 'un', false);
%                                         
%     scaled_threat_intensity = cellfun(@(x) vertcat(x{:}), scaled_threat_intensity, 'un', false);
%     
%     merged_finalsale_data.scaled_threat_intensity = vertcat(scaled_threat_intensity{:});
%     
% end

% 
% function merged_aggregated_sector_scale = merge_aggregated_sector_scale_data_test(consumption_footprint_data, finalsale_data, industry_characteristics, analyse_footprint_params)
%     
%     finalsale_link_indexes = arrayfun(@(x) find(consumption_footprint_data.aggregated_sector_scale.finalsale_sector_list == x), finalsale_data.aggregated_sector_scale.aggregated_paths, 'un', false);
% 
%     non_empties = cellfun('length', finalsale_link_indexes) > 0;
%     
%     finalsale_link_indexes = finalsale_link_indexes(non_empties);
%     finalsale_threat_intensities_to_use = finalsale_data.aggregated_sector_scale.aggregated_path_vals(non_empties);
% 
%     consumption_country_names_per_finalsale_commodity = cellfun(@(x) consumption_footprint_data.aggregated_sector_scale.consumption_country_list(x), finalsale_link_indexes, 'un', false);
% 
%     expanded_finalsale_commodity_names = cellfun(@(x) industry_characteristics.commodity_classification_list(consumption_footprint_data.aggregated_sector_scale.finalsale_sector_list(x)), ...
%                                                                         finalsale_link_indexes, 'un', false);
% 
%     expanded_finalsale_country_names = cellfun(@(x) industry_characteristics.country_names_list(consumption_footprint_data.aggregated_sector_scale.finalsale_sector_list(x)), ...
%                                                                         finalsale_link_indexes, 'un', false);    
%                                                                     
%     merged_aggregated_sector_scale.scaled_threat_proportion = cellfun(@(x, y) consumption_footprint_data.aggregated_sector_scale.threat_intensities(x) ./ (sum(consumption_footprint_data.aggregated_sector_scale.threat_intensities(x))), ...
%                                                                        finalsale_link_indexes, 'un', false);
% 
%     merged_aggregated_sector_scale.scaled_threat_intensity = cellfun(@(x, y) x.*y, merged_aggregated_sector_scale.scaled_threat_proportion, num2cell(finalsale_threat_intensities_to_use), 'un', false);
%     
%     merged_aggregated_sector_scale.scaled_threat_intensity = vertcat(merged_aggregated_sector_scale.scaled_threat_intensity{:});
%     
%     expanded_net_threat_intensities = cellfun(@(x, y) num2cell(repmat(finalsale_threat_intensities_to_use(y), [size(x, 1) 1])), ...
%                                                                             finalsale_link_indexes, num2cell(1:length(finalsale_link_indexes))', 'un', false);
%                                                                         
%     expanded_net_threat_intensities = vertcat(merged_aggregated_sector_scale.expanded_net_threat_intensities{:});
%     
%     expanded_species_counts = cell2struct(cell(length(analyse_footprint_params.groups_to_count), 1), analyse_footprint_params.groups_to_count);                                       
%     
%     for group_ind = 1:length(analyse_footprint_params.groups_to_count)
%         expanded_species_counts.(analyse_footprint_params.groups_to_count{group_ind}) = cellfun(@(x, y) num2cell(repmat(finalsale_data.aggregated_sector_scale.species_counts.(analyse_footprint_params.groups_to_count{group_ind})(y), [size(x, 1) 1])), ...
%                                                                                             finalsale_link_indexes, num2cell(1:length(finalsale_link_indexes))', 'un', false); 
%     end  
%     
%     merged_aggregated_sector_scale.data_block = [vertcat(consumption_country_names_per_finalsale_commodity{:})...
%                                                  vertcat(expanded_finalsale_country_names{:})...
%                                                  vertcat(expanded_finalsale_commodity_names{:})...
%                                                  vertcat(expanded_species_counts.(analyse_footprint_params.groups_to_count{1}){:}) ...
%                                                  vertcat(expanded_species_counts.(analyse_footprint_params.groups_to_count{2}){:}) ...
%                                                  num2cell(merged_aggregated_sector_scale.scaled_threat_intensity), ...
%                                                  num2cell(vertcat(merged_aggregated_sector_scale.scaled_threat_proportion{:})) ...
%                                                  vertcat(expanded_net_threat_intensities{:})];
%  
%     [~, indexes_sorted_by_threat_intensity] = sort(merged_aggregated_sector_scale.scaled_threat_intensity, 'descend');
%     
%     merged_aggregated_sector_scale.table = cell2table([num2cell(1:length(merged_aggregated_sector_scale.scaled_threat_intensity))'...
%                                                                         merged_aggregated_sector_scale.data_block(indexes_sorted_by_threat_intensity, :)], ...
%                                                              'VariableNames', [{'global_rank', 'consumption_country', 'finalsale_country', 'finalsale_industry'}, analyse_footprint_params.groups_to_count, ...
%                                                                                {'attributed_threat_intensity', 'finalsale_threat_intensity_proportion', 'net_threat_intensity'}]);
%                                                               
% end
% 


% function processed_finalsale_data = process_footprint_at_finalsale_level(industry_characteristics, analyse_footprint_params, footprint_filename_prefix, species_characteristics, species_to_use, species_group_weight)
%      
%     
%    % load([footprint_filename_prefix, 'finalsale_', num2str(species_to_use(species_ind)), '.mat'])
%    
%     finalsale_footprints = load([footprint_filename_prefix, 'finalsale_', num2str(species_to_use(species_ind)), '.mat']);
%     
%     if analyse_footprint_params.weight_species == true
%         finalsale_footprints.vals = species_group_weight*finalsale_footprints.vals/numel(species_to_use);
%     end
%     
%     finalsale_footprints = select_footprint_subset(finalsale_footprints, ...
%                                                    analyse_footprint_params, ...
%                                                    species_to_use, ...
%                                                    industry_characteristics);                                                                                       
% 
%     if ~isempty(finalsale_footprints.species)                                                
%         
%         processed_finalsale_data = aggregate_footprint_at_finalsale_routines(finalsale_footprints, analyse_footprint_params, industry_characteristics);
%                      
%         processed_finalsale_data.aggregated_sector_scale = structfun(@(x) x, processed_finalsale_data.aggregated_sector_scale, 'un', false);
%     
%         display([sum(processed_finalsale_data.country_scale.finalsale.aggregated_path_vals)])
%         processed_finalsale_data.aggregated_sector_scale.species_counts = count_species_groups(processed_finalsale_data.aggregated_sector_scale.grouped_aggregates, species_characteristics, analyse_footprint_params.groups_to_count);
%     else 
%         processed_finalsale_data = struct();
%     end
%     
% end
% 
% 
% function consumption_finalsale_matches = match_consumption_to_finalsale(consumption_footprint_data, processed_finalsale_data, country_indexes_to_use)
%     
%     consumption_finalsale_matches = struct();
%                                                                   
%     [consumption_finalsale_matches.consumption_commodity, consumption_index, ~] = arrayfun(@(x) intersect(grouped_finalsale_data{x}.aggregated_sector_scale.aggregated_paths, ...
%                                                                                                           processed_finalsale_data.aggregated_sector_scale.aggregated_paths), (1:numel(country_indexes_to_use))', 'un', false);
%     
%     consumption_finalsale_matches.consumption_threat_intensities = cellfun(@(x, y) grouped_finalsale_data{x}.aggregated_sector_scale.aggregated_path_vals(y), ...
%                                                                                     num2cell(1:numel(country_indexes_to_use))', consumption_index, 'un', false);
%                             
%     consumption_finalsale_matches.consumption_country_index = cellfun(@(x, y) repmat(x, [numel(y) 1]), num2cell(country_indexes_to_use), consumption_finalsale_matches.consumption_commodity, 'un', false);
% 
%     consumption_finalsale_matches = structfun(@(x) vertcat(x{:}), consumption_finalsale_matches, 'un', false);
% 
% end
% 
% 
% 
% 
% 
% function merged_finalsale_data = merge_consumption_data_with_finalsale(consumption_finalsale_matches, finalsale_data, analyse_footprint_params, industry_characteristics)
% 
%     consumption_country_grouped_by_commodity = aggregate_and_sort_paths(analyse_footprint_params.sort_data, ...
%                                                                         consumption_finalsale_matches.consumption_threat_intensities, ...
%                                                                         consumption_finalsale_matches.consumption_commodity, ...
%                                                                         consumption_finalsale_matches.consumption_country_index);    
%                                                                                                 
%     [~, ~, matched_consumption_finalsale_commodity_indexes] = intersect(consumption_country_grouped_by_commodity.aggregated_paths, finalsale_data.aggregated_sector_scale.aggregated_paths, 'stable');
%         
%     merged_finalsale_data.scaled_threat_intensity_blocks = build_scaled_threat_intensities(consumption_country_grouped_by_commodity, ...
%                                                                                             finalsale_data,matched_consumption_finalsale_commodity_indexes,...
%                                                                                             finalsale_data.aggregated_sector_scale.aggregated_path_vals(matched_consumption_finalsale_commodity_indexes));
%      
%     merged_finalsale_data.consumption_finalsale_blocks = build_consumption_finalsale_production_blocks(consumption_country_grouped_by_commodity.grouped_aggregates, ...
%                                                                                                        finalsale_data.aggregated_sector_scale.aggregated_paths(matched_consumption_finalsale_commodity_indexes), ...
%                                                                                                        finalsale_data.aggregated_sector_scale.grouped_sectors(matched_consumption_finalsale_commodity_indexes));
% 
%                                                                                        
%     merged_finalsale_data.international_indexes = sum( abs(diff([merged_finalsale_data.consumption_finalsale_blocks.consumption_country, ...
%                                                                 industry_characteristics.country_index_list(merged_finalsale_data.consumption_finalsale_blocks.finalsale), ...
%                                                                 industry_characteristics.country_index_list(merged_finalsale_data.consumption_finalsale_blocks.production)], 1, 2)), 2) > 0;  
%                                                             
% end
% 



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD TABLE BUILD ROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% function T = build_aggregated_table(rank_type, groups_to_count, industry_characteristics, object_to_aggregate_over, net_threat_aggregate, species_counts)
%     
%     if strcmp(rank_type, 'by_industry')
%         aggregated_outputs = cell(length(object_to_aggregate_over), 2);
%         aggregated_outputs(:, 1) = industry_characteristics.country_names_list( object_to_aggregate_over);
%         aggregated_outputs(:, 2) = industry_characteristics.commodity_classification_list(object_to_aggregate_over);
%         table_names = [{'Consumption_Country', 'Consumption_Industry'}, groups_to_count, {'Aggregated_Threats'}];
%         
%     else 
%         aggregated_outputs = cell(length(object_to_aggregate_over), 1);
%         aggregated_outputs(:, 1) = industry_characteristics.country_index_map( object_to_aggregate_over, 1);
%         table_names = [{'Consumption_Country'}, groups_to_count, {'Aggregated_Threats'}];
% 
%     end
% 
%     aggregated_outputs = [aggregated_outputs num2cell([species_counts{:}])  num2cell(round(net_threat_aggregate, 1))];
%     
%     aggregated_outputs = append_zero_block(aggregated_outputs, industry_characteristics, object_to_aggregate_over);
%         
%     T = cell2table(aggregated_outputs, 'VariableNames', table_names);
%     
% end
% 
% 
% 
% function expanded_mrio_data = write_expanded_ranked_outputs(ranked_mrio_industries, analyse_footprint_params, industry_characteristics, species_characteristics, grouped_mrio_industry_paths, industry_grouped_mrio_species_aggregates, grouped_mrio_industry_vals)
% 
%         ranked_mrio_industries = arrayfun(@(x) ranked_mrio_industries(x, :), (1:length(ranked_mrio_industries))', 'un', 0);
%         ranked_mrio_industries = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), size(x, 2))), ranked_mrio_industries, grouped_mrio_industry_paths, 'un', 0);
%         ranked_mrio_industries = vertcat(ranked_mrio_industries{:});
%     
%         expanded_mrio_data = cell(1, 5);
%         expanded_mrio_data{1} = cellfun(@(x) industry_characteristics.country_names_list(x), grouped_mrio_industry_paths, 'un', false); 
%         expanded_mrio_data{2} = cellfun(@(x) industry_characteristics.commodity_classification_list(x), grouped_mrio_industry_paths, 'un', false); 
%         expanded_mrio_data{3} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'ANIMALIA'))), x,'un', 0), industry_grouped_mrio_species_aggregates,'un', 0);
%         expanded_mrio_data{4} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'PLANTAE'))), x,'un', 0), industry_grouped_mrio_species_aggregates,'un', 0);
%         %expanded_mrio_data{5} = cellfun(@(x) cellfun(@(y) length(y), x,'un',0), industry_grouped_mrio_species_aggregates,'un', 0);
%         expanded_mrio_data{5} = cellfun(@(x) num2cell(round(x, 1)), grouped_mrio_industry_vals, 'un', 0);
%         
%         %expanded_mrio_data{7} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'FUNGI'))), x,'un', 0), industry_grouped_mrio_species_aggregates,'un', 0);
%         %expanded_mrio_data{8} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'CHROMISTA'))), x,'un', 0), industry_grouped_mrio_species_aggregates,'un', 0);
%     
%         expanded_mrio_data = cellfun(@(x) vertcat(x{:}), expanded_mrio_data, 'un', 0);
%         expanded_mrio_data = [expanded_mrio_data{:}];
%         expanded_mrio_data = [ranked_mrio_industries expanded_mrio_data];
%     
%         T = cell2table(expanded_mrio_data, 'VariableNames', {'Consumption_Country', 'Consumption_Industry', 'Total_Threatened_Animalia', ...
%                                                             'Total_Threatened_Plantae',  'Aggregated_Threats', 'Production_Country', 'Production_Industry', ...
%                                                             'Threatened_Animalia', 'Threatened_Plantae','Aggregated_Threats_Per_Tradepath' });
%                                                         
%         writetable(T, strcat(analyse_footprint_params.output_folder, analyse_footprint_params.country_of_interest, '_', analyse_footprint_params.impact_assessment_level, '_', analyse_footprint_params.finalsale_scale, '_finalsale_scale_expanded.txt'), 'Delimiter', 'tab')
%     
% end


% function T = write_sector_to_sector_scale_table(industry_characteristics, species_characteristics, ranked_industry_paths_per_species, ranked_aggregated_threats_per_species, grouped_mrio_species_aggregates)
%        
%     aggregated_outputs = cell(length(ranked_industry_paths_per_species), 4);
%     aggregated_outputs(:, 1) = industry_characteristics.country_names_list( ranked_industry_paths_per_species(:, 1)); 
%     aggregated_outputs(:, 2) = industry_characteristics.commodity_classification_list(ranked_industry_paths_per_species(:, 1));
%     aggregated_outputs(:, 3) = industry_characteristics.country_names_list( ranked_industry_paths_per_species(:, 2)); 
%     aggregated_outputs(:, 4) = industry_characteristics.commodity_classification_list(ranked_industry_paths_per_species(:, 2));
%     
%     decomposed_aggregate_species_characteristics = cell(length(grouped_mrio_species_aggregates), 2);
% 
%     decomposed_aggregate_species_characteristics(:, 1) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'ANIMALIA')))), grouped_mrio_species_aggregates));
%     decomposed_aggregate_species_characteristics(:, 2) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'PLANTAE')))), grouped_mrio_species_aggregates));
%     
%     %decomposed_aggregate_species_characteristics(:, 3) = num2cell(cellfun(@(x) length(unique(x)), grouped_mrio_species_aggregates));
%     %decomposed_aggregate_species_characteristics(:, 4) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'FUNGI')))), grouped_mrio_species_aggregates));
%     %decomposed_aggregate_species_characteristics(:, 5) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'CHROMISTA')))), grouped_mrio_species_aggregates));
%     
%     aggregated_outputs = [aggregated_outputs decomposed_aggregate_species_characteristics num2cell(round(ranked_aggregated_threats_per_species, 1))];
%     
%     table_names = {'Consumption_Country', 'Consumption_Industry', 'Production_Country', 'Production_Industry', 'Threatened_Animalia',  'Threatened_Plantae', 'Aggregated_Threats'};
%        
%     T = cell2table(aggregated_outputs, 'VariableNames', table_names);
%     
% end
% 
% 
% function ranked_countries = append_zero_block(ranked_countries, industry_characteristics, country_indexes_used)
%     consumption_countries_to_append = setdiff(cell2mat(industry_characteristics.country_index_map(:, 2)), country_indexes_used);
%         
%     cellblock_to_append = cell(length(consumption_countries_to_append), size(ranked_countries, 2));
%     cellblock_to_append(:, 1) = industry_characteristics.country_index_map(consumption_countries_to_append, 1);
%         
%     cellblock_to_append(:, 2:end) = num2cell(zeros(size(cellblock_to_append, 1), (size(cellblock_to_append, 2) - 1)));
%     
%     ranked_countries = vertcat(ranked_countries, cellblock_to_append);
% end


