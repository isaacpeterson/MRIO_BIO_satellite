function balanced_footprint = build_balanced_footprint(consumption_footprint_characteristics, footprint_input_objects, industry_inputs)

    balanced_footprint = struct();
    balanced_footprint.sector_to_sector_scale = aggregate_grouped_data(consumption_footprint_characteristics, footprint_input_objects, ...
                                                                                                        industry_inputs.industry_characteristics, 'sector_to_sector_scale');
                                                                                                    
    balanced_footprint.aggregated_sector_scale = aggregate_grouped_data(consumption_footprint_characteristics, footprint_input_objects, ...
                                                                                                        industry_inputs.industry_characteristics, 'aggregated_sector_scale');
    balanced_footprint.country_scale = aggregate_grouped_data(consumption_footprint_characteristics, footprint_input_objects, ...
                                                                                                        industry_inputs.industry_characteristics, 'country_scale'); 
                                                                                                    
end
function balanced_data = aggregate_grouped_data(consumption_footprint_characteristics, footprint_input_objects, industry_characteristics, trade_scale)
    
    countries_to_group = structfun(@(x) x.(trade_scale).consumption_country, consumption_footprint_characteristics, 'un', false);
    
    if strcmp(trade_scale, 'country_scale')
        paths_to_group = struct2cell(countries_to_group);
    else 
        paths_to_group = structfun(@(x) x.(trade_scale).paths, consumption_footprint_characteristics, 'un', false);
        paths_to_group = cellfun(@(x) horzcat(countries_to_group.(x), paths_to_group.(x)), footprint_input_objects.species_group_names, 'un', false);
    end
    
    unique_characteristics = unique(vertcat(paths_to_group{:}), 'rows');
    
    vals_to_aggregate = structfun(@(x) x.(trade_scale).threat_intensities, consumption_footprint_characteristics, 'un', false);
    vals_to_aggregate = struct2cell(vals_to_aggregate);
    
    sorted_data = cellfun(@(x) balance_datablock(vals_to_aggregate{x}, paths_to_group{x}, unique_characteristics), num2cell(1:4)', 'un', false);
            
    scaled_data = cellfun(@(x, y) x./numel(y), sorted_data, footprint_input_objects.satellite_collapse_groups, 'un', false);
    
    aggregated_data = sum(cat(3, sorted_data{:}), 3);
    scaled_aggregated_data = sum(cat(3, scaled_data{:}), 3);
    
    path_block = industry_characteristics.unique_countries(unique_characteristics(:, 1));
    path_names = {'Consumption_Country'};
    
    if ~strcmp(trade_scale, 'country_scale')
        
        if strcmp(trade_scale, 'sector_to_sector_scale')
        	path_names = [path_names {'finalsale_country', 'finalsale_industry', 'production_country', 'production_industry'}];
            industry_block = [industry_characteristics.country_names_list(unique_characteristics(:, 2)), ...
                     industry_characteristics.commodity_classification_list(unique_characteristics(:, 2))...
                     industry_characteristics.country_names_list(unique_characteristics(:, 3)), ...
                     industry_characteristics.commodity_classification_list(unique_characteristics(:, 3))];
        
        elseif strcmp(trade_scale, 'aggregated_sector_scale')
            path_names = [path_names {'finalsale_country', 'finalsale_industry'}];
            industry_block = [industry_characteristics.country_names_list(unique_characteristics(:, 2)), ...
                     industry_characteristics.commodity_classification_list(unique_characteristics(:, 2))];
                 
        end
        
        path_block = [path_block industry_block]; 
    end
    
    [~, sorted_inds] = sort(scaled_aggregated_data, 'descend');
    data_block = [path_block num2cell([sorted_data{:} aggregated_data scaled_aggregated_data])];
    balanced_data = cell2table(data_block(sorted_inds, :), 'VariableNames', [path_names, footprint_input_objects.species_group_names{:}, {'total', 'balanced'}]);
    
    
    %     norm_const = sum(sorted_data(:, strcmp(current_group.datablock_names, 'total')));
%     
%     groups_to_balance = ~ismember(current_group.datablock_names, species_group_names);
%     
%     sorted_data(:, groups_to_balance) = current_weight * sorted_data(:, groups_to_balance)/norm_const;
end


function sorted_data = balance_datablock(block_to_aggregate, current_characteristic_set, unique_characteristics)
    
    sorted_data = zeros(size(unique_characteristics, 1), size(block_to_aggregate, 2));
    [~, template_indexes, sorted_data_indexes] = intersect(unique_characteristics, current_characteristic_set, 'rows', 'stable');
    
    sorted_data(template_indexes, :) = block_to_aggregate(sorted_data_indexes, :);
    
end

