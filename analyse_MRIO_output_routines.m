function [outputs] = analyse_MRIO_output_routines(analyse_MRIO_params, IUCN_data_object, MRIO_threat_tensor, species_characteristics)
    
    tic
    
    outputs = struct();
    
    [MRIO_species_industry_identifiers, MRIO_species_threat_proportions] = select_MRIO_subset(analyse_MRIO_params, MRIO_threat_tensor, species_characteristics, ...
                                                                                            IUCN_data_object.industry_characteristics, IUCN_data_object.IUCN_threat_taxons, IUCN_data_object.IUCN_threat_status);
                                                                                                      
    species_level = aggregate_and_sort_paths(analyse_MRIO_params.sort_data, analyse_MRIO_params.groups_to_count, MRIO_species_threat_proportions, MRIO_species_industry_identifiers(:, 1:2), MRIO_species_industry_identifiers(:, 3), species_characteristics, false);
   
    
    if ~isempty(fieldnames(species_level))
        outputs.species_level = species_level;
        outputs.species_level.table = write_species_level_table(IUCN_data_object.industry_characteristics, species_characteristics, outputs.species_level.aggregated_paths, outputs.species_level.aggregated_path_vals, outputs.species_level.grouped_aggregates);

        outputs.industry_level = aggregate_and_sort_paths(analyse_MRIO_params.sort_data, analyse_MRIO_params.groups_to_count, outputs.species_level.aggregated_path_vals, outputs.species_level.aggregated_paths(:, 1), outputs.species_level.aggregated_paths(:, 2), species_characteristics, false);
        outputs.industry_level.aggregates = aggregate_and_sort_paths(analyse_MRIO_params.sort_data, analyse_MRIO_params.groups_to_count, outputs.species_level.aggregated_path_vals, outputs.species_level.aggregated_paths(:, 1), outputs.species_level.grouped_aggregates, species_characteristics, true);
   
        outputs.industry_level.table = build_aggregated_table('by_industry', analyse_MRIO_params.groups_to_count, IUCN_data_object.industry_characteristics, outputs.industry_level.aggregated_paths, ...
                                                            outputs.industry_level.aggregates.aggregated_path_vals, outputs.industry_level.aggregates.species_counts);
        
        outputs.country_level.production = aggregate_and_sort_paths(analyse_MRIO_params.sort_data, analyse_MRIO_params.groups_to_count, outputs.species_level.aggregated_path_vals, ...
                                                                    IUCN_data_object.industry_characteristics.country_index_list(outputs.species_level.aggregated_paths(:, 2)), outputs.species_level.grouped_aggregates, species_characteristics, true);
  
        outputs.country_level.production.table = build_aggregated_table('by_country', analyse_MRIO_params.groups_to_count, IUCN_data_object.industry_characteristics, outputs.country_level.production.aggregated_paths, ... 
                                                         outputs.country_level.production.aggregated_path_vals, outputs.country_level.production.species_counts);
   
        outputs.country_level.finalsale = aggregate_and_sort_paths(analyse_MRIO_params.sort_data, analyse_MRIO_params.groups_to_count, outputs.species_level.aggregated_path_vals, IUCN_data_object.industry_characteristics.country_index_list(outputs.species_level.aggregated_paths(:, 1)), outputs.species_level.grouped_aggregates, species_characteristics, true);
        outputs.country_level.finalsale.table = build_aggregated_table('by_country', analyse_MRIO_params.groups_to_count, IUCN_data_object.industry_characteristics, outputs.country_level.finalsale.aggregated_paths, ... 
                                                                    outputs.country_level.finalsale.aggregated_path_vals, outputs.country_level.finalsale.species_counts);        
    end
%     if analyse_MRIO_params.write_industry_ranks == true                                                           
%         filename = strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_species_level_', analyse_MRIO_params.assessment_scale, '_assessment_scale.txt');
%         writetable(outputs.species_level.table, filename, 'Delimiter', 'tab')
%     end
%     if analyse_MRIO_params.write_production_country_ranks == true                                                      
%         filename = strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_', 'production', '_', analyse_MRIO_params.assessment_scale, '_assessment_scale.txt');
%         writetable(outputs.industry_level.table, filename, 'Delimiter', 'tab')                                                  
%     end
%                                                            
%     if analyse_MRIO_params.write_finalsale_country_ranks == true
%         
%     end
%     
%     if analyse_MRIO_params.write_industry_table == true
%         filename = strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_', file_identifier, '_', analyse_MRIO_params.assessment_scale, '_industry_assessment_scale.txt');
%         writetable(outputs.industry_level.table, filename, 'Delimiter', 'tab')
% %         if analyse_MRIO_params.write_expanded_table == true
% %             outputs.ranked_expanded_paths = write_expanded_ranked_outputs(outputs.ranked_aggregated_paths, analyse_MRIO_params, IUCN_data_object.industry_characteristics, species_characteristics, outputs.industry_level.grouped_aggregates, outputs.industry_level.aggregates.grouped_aggregates, outputs.industry_level.grouped_path_vals);
% %         end
%     end
end


function T = write_species_level_table(industry_characteristics, species_characteristics, ranked_industry_paths_per_species, ranked_aggregated_threats_per_species, grouped_MRIO_species_aggregates)
       
    aggregated_outputs = cell(length(ranked_industry_paths_per_species), 4);
    aggregated_outputs(:, 1) = industry_characteristics.country_names_list( ranked_industry_paths_per_species(:, 1)); 
    aggregated_outputs(:, 2) = industry_characteristics.commodity_classification_list(ranked_industry_paths_per_species(:, 1));
    aggregated_outputs(:, 3) = industry_characteristics.country_names_list( ranked_industry_paths_per_species(:, 2)); 
    aggregated_outputs(:, 4) = industry_characteristics.commodity_classification_list(ranked_industry_paths_per_species(:, 2));
    
    decomposed_aggregate_species_characteristics = cell(length(grouped_MRIO_species_aggregates), 2);

    decomposed_aggregate_species_characteristics(:, 1) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'ANIMALIA')))), grouped_MRIO_species_aggregates));
    decomposed_aggregate_species_characteristics(:, 2) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'PLANTAE')))), grouped_MRIO_species_aggregates));
    
    %decomposed_aggregate_species_characteristics(:, 3) = num2cell(cellfun(@(x) length(unique(x)), grouped_MRIO_species_aggregates));
    %decomposed_aggregate_species_characteristics(:, 4) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'FUNGI')))), grouped_MRIO_species_aggregates));
    %decomposed_aggregate_species_characteristics(:, 5) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'CHROMISTA')))), grouped_MRIO_species_aggregates));
    
    
    aggregated_outputs = [aggregated_outputs decomposed_aggregate_species_characteristics num2cell(round(ranked_aggregated_threats_per_species, 1))];
    
    table_names = {'Consumption_Country', 'Consumption_Industry', 'Production_Country', 'Production_Industry', 'Threatened_Animalia',  'Threatened_Plantae', 'Aggregated_Threats'};
       
    T = cell2table(aggregated_outputs, 'VariableNames', table_names);
    
end



function ranked_countries = append_zero_block(ranked_countries, industry_characteristics, country_indexes_used)
    consumption_countries_to_append = setdiff(cell2mat(industry_characteristics.country_index_map(:, 2)), country_indexes_used);
        
    cellblock_to_append = cell(length(consumption_countries_to_append), size(ranked_countries, 2));
    cellblock_to_append(:, 1) = industry_characteristics.country_index_map(consumption_countries_to_append, 1);
        
    cellblock_to_append(:, 2:end) = num2cell(zeros(size(cellblock_to_append, 1), (size(cellblock_to_append, 2) - 1)));
    
    ranked_countries = vertcat(ranked_countries, cellblock_to_append);
end

function species_group_counts = count_species_groups(industry_grouped_species_aggregates, species_characteristics, groups_to_count)
    
    species_group_counts = cell(1, length(groups_to_count));
    
    for current_counter = 1:length(groups_to_count)
        %species_group_counts{current_counter} = cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), groups_to_count{current_counter})), x, 'un', 0) ))), industry_grouped_species_aggregates);
%         species_to_use = cellfun(@(x) vertcat(x{:}), industry_grouped_species_aggregates, 'un', false);
%         species_group_counts{current_counter} = cell2mat(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), groups_to_count{current_counter})))), species_to_use, 'un', false));
        species_to_use = cellfun(@(x) unique(vertcat(x{:})), industry_grouped_species_aggregates, 'un', false);
        species_group_counts{current_counter} = cell2mat(cellfun(@(x) length(find(strcmp(species_characteristics.species_kingdom(x), groups_to_count{current_counter}))), species_to_use, 'un', false));
    end
    
end

function T = build_aggregated_table(rank_type, groups_to_count, industry_characteristics, object_to_aggregate_over, net_threat_aggregate, species_counts)
    
    if strcmp(rank_type, 'by_industry')
        aggregated_outputs = cell(length(object_to_aggregate_over), 2);
        aggregated_outputs(:, 1) = industry_characteristics.country_names_list( object_to_aggregate_over);
        aggregated_outputs(:, 2) = industry_characteristics.commodity_classification_list(object_to_aggregate_over);
        table_names = [{'Consumption_Country', 'Consumption_Industry'}, groups_to_count, {'Aggregated_Threats'}];
        
    else 
        aggregated_outputs = cell(length(object_to_aggregate_over), 1);
        aggregated_outputs(:, 1) = industry_characteristics.country_index_map( object_to_aggregate_over, 1);
        table_names = [{'Consumption_Country'}, groups_to_count, {'Aggregated_Threats'}];

    end

    aggregated_outputs = [aggregated_outputs num2cell([species_counts{:}])  num2cell(round(net_threat_aggregate, 1))];
    
    aggregated_outputs = append_zero_block(aggregated_outputs, industry_characteristics, object_to_aggregate_over);
        
    T = cell2table(aggregated_outputs, 'VariableNames', table_names);
end



function expanded_MRIO_data = write_expanded_ranked_outputs(ranked_MRIO_industries, analyse_MRIO_params, industry_characteristics, species_characteristics, grouped_MRIO_industry_paths, industry_grouped_MRIO_species_aggregates, grouped_MRIO_industry_vals)

        ranked_MRIO_industries = arrayfun(@(x) ranked_MRIO_industries(x, :), (1:length(ranked_MRIO_industries))', 'un', 0);
        ranked_MRIO_industries = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), size(x, 2))), ranked_MRIO_industries, grouped_MRIO_industry_paths, 'un', 0);
        ranked_MRIO_industries = vertcat(ranked_MRIO_industries{:});
    
        expanded_MRIO_data = cell(1, 5);
        expanded_MRIO_data{1} = cellfun(@(x) industry_characteristics.country_names_list(x), grouped_MRIO_industry_paths, 'un', false); 
        expanded_MRIO_data{2} = cellfun(@(x) industry_characteristics.commodity_classification_list(x), grouped_MRIO_industry_paths, 'un', false); 
        expanded_MRIO_data{3} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'ANIMALIA'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
        expanded_MRIO_data{4} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'PLANTAE'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
        %expanded_MRIO_data{5} = cellfun(@(x) cellfun(@(y) length(y), x,'un',0), industry_grouped_MRIO_species_aggregates,'un', 0);
        expanded_MRIO_data{5} = cellfun(@(x) num2cell(round(x, 1)), grouped_MRIO_industry_vals, 'un', 0);
        
        %expanded_MRIO_data{7} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'FUNGI'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
        %expanded_MRIO_data{8} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'CHROMISTA'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
    
        expanded_MRIO_data = cellfun(@(x) vertcat(x{:}), expanded_MRIO_data, 'un', 0);
        expanded_MRIO_data = [expanded_MRIO_data{:}];
        expanded_MRIO_data = [ranked_MRIO_industries expanded_MRIO_data];
    
        T = cell2table(expanded_MRIO_data, 'VariableNames', {'Consumption_Country', 'Consumption_Industry', 'Total_Threatened_Animalia', ...
                                                            'Total_Threatened_Plantae',  'Aggregated_Threats', 'Production_Country', 'Production_Industry', ...
                                                            'Threatened_Animalia', 'Threatened_Plantae','Aggregated_Threats_Per_Tradepath' });
        writetable(T, strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_', analyse_MRIO_params.industry_assessment_type, '_', analyse_MRIO_params.assessment_scale, '_assessment_scale_expanded.txt'), 'Delimiter', 'tab')
    
end



%(analyse_MRIO_params.sort_data, analyse_MRIO_params.groups_to_count, outputs.species_level.aggregated_path_vals, outputs.species_level.aggregated_paths(:, 1), outputs.species_level.grouped_aggregates, species_characteristics, true);
   

function [outputs] = aggregate_and_sort_paths(sort_data, groups_to_count, vals_to_aggregate, paths_to_aggregate, objects_to_aggregate, species_characteristics, calc_species_num)
    
    outputs = struct();
    [aggregated_paths, ~, path_indexes] = unique(paths_to_aggregate, 'rows', 'stable');
    if ~isempty(path_indexes)
        
        outputs.aggregated_paths = aggregated_paths;
        outputs.aggregated_path_vals = accumarray(path_indexes, vals_to_aggregate);
    

        grouped_path_indexes = accumarray(path_indexes, find(path_indexes), [], @(rows){rows}); 
         if sort_data == true
            [outputs.aggregated_path_vals, sorted_indexes] = sort(outputs.aggregated_path_vals, 'descend');
            outputs.aggregated_paths = outputs.aggregated_paths(sorted_indexes, :);
            grouped_path_indexes = grouped_path_indexes(sorted_indexes);
            grouped_path_indexes = cellfun(@(x) x(sort_indexes(vals_to_aggregate(x))), grouped_path_indexes, 'un', false);
         end
    
        outputs.grouped_path_vals = cellfun(@(x) vals_to_aggregate(x), grouped_path_indexes, 'un', false);
        outputs.grouped_aggregates = cellfun(@(x) objects_to_aggregate(x, :), grouped_path_indexes, 'un', false);
    
        if calc_species_num == true     
            outputs.species_counts = count_species_groups(outputs.grouped_aggregates, species_characteristics, groups_to_count);
        end

    end
   
    %outputs.species_level.aggregated_paths(:, 1), outputs.species_level.aggregated_paths(:, 2)
end

function [sorted_inds] = sort_indexes(arr)
    [~, sorted_inds] = sort(arr, 'descend');
end

% 
% function [domestic_threat_paths, international_threat_paths, global_threat_paths] = build_MRIO_threat_rankings(MRIO_consumption_production, species_characteristics, analyse_MRIO_params, ...
%                                                                                                             MRIO_species_threat_proportions, MRIO_species_industry_identifiers, IUCN_species_taxons)
%      
%     MRIO_threat_set = build_MRIO_threat_set(species_characteristics.species_taxons(MRIO_species_industry_identifiers(:, 3)), IUCN_species_taxons);
%     
%     MRIO_consumption_production{5} = species_characteristics.species_names(MRIO_species_industry_identifiers(:, 3));
%     MRIO_consumption_production{6} = num2cell(MRIO_threat_set);
%     MRIO_consumption_production{7} = num2cell(MRIO_species_threat_proportions);
%     global_threat_paths = [MRIO_consumption_production{:}];
%     [~, sorted_indexes] = sort(MRIO_species_threat_proportions, 'descend');
%     global_threat_paths = global_threat_paths(sorted_indexes, :);
%      
%     if strcmp(analyse_MRIO_params.industry_assessment_type, 'production_based')
%         domestic_threat_indexes = strcmp(MRIO_consumption_production{1}, analyse_MRIO_params.country_of_interest);
%     elseif strcmp(analyse_MRIO_params.industry_assessment_type, 'finalsale_based')
%         domestic_threat_indexes = strcmp(MRIO_consumption_production{3}, analyse_MRIO_params.country_of_interest);
%     end
% 
%     domestic_threat_paths = global_threat_paths(domestic_threat_indexes, :);
%     international_threat_paths = global_threat_paths(~domestic_threat_indexes, :);
%     
% end
%     
% function save_MRIO_threat_rankings(domestic_threat_paths, international_threat_paths, global_threat_paths)
%     if strcmp(analyse_MRIO_params.industry_assessment_type, 'production_based')
%         
%         file_prefix = strcat(analyse_MRIO_params.output_folder, 'threat_proportion_rankings_consumption_type_');
%         file_suffix = strcat('_production_country_', analyse_MRIO_params.country_of_interest, '.mat');
%         int_filename = strcat(file_prefix, 'international', file_suffix);
%         dom_filename = strcat(file_prefix, 'domestic', file_suffix);
%         global_filename = strcat(file_prefix, 'global', file_suffix);
%     elseif strcmp(analyse_MRIO_params.industry_assessment_type, 'finalsale_based')
%         
%         file_prefix = strcat(analyse_MRIO_params.output_folder, 'threat_proportion_rankings_consumption_country_', analyse_MRIO_params.country_of_interest, '_production_type_');
%         file_suffix =  strcat('.mat');
%         int_filename = strcat(file_prefix, 'international', file_suffix);
%         dom_filename = strcat(file_prefix, 'domestic' , file_suffix);
%         global_filename = strcat(file_prefix, 'global', file_suffix);
%     end
%              
%     save(int_filename, 'international_threat_paths', '-v7.3')
%     save(dom_filename, 'domestic_threat_paths', '-v7.3')
%     save(global_filename, 'global_threat_paths', '-v7.3')
% end

function MRIO_threat_set = build_MRIO_threat_set(MRIO_species_taxons, IUCN_species_taxons) 
    
    unique_IUCN_species_taxons = unique(IUCN_species_taxons);
    [species_threat_num] = histc(IUCN_species_taxons, unique_IUCN_species_taxons);

    [unique_MRIO_species, ~, threat_indexes] = intersect(MRIO_species_taxons, unique_IUCN_species_taxons, 'stable');

    MRIO_threat_set = zeros(length(MRIO_species_taxons), 1);

    for current_ind = 1:length(threat_indexes)
        MRIO_threat_set(MRIO_species_taxons == unique_MRIO_species(current_ind)) = species_threat_num(threat_indexes(current_ind));  
    end
    
end

function industry_characteristics = build_industry_characteristics(analyse_MRIO_params)
    
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

end



% 
% 
% function [MRIO_species_industry_identifiers, MRIO_species_threat_proportions] = select_MRIO_subset(analyse_MRIO_params, MRIO_threat_tensor, species_characteristics, ...
%                                                                                                 industry_characteristics, IUCN_threat_taxons, IUCN_threat_status)
%     
%     if (strcmp(analyse_MRIO_params.status_levels_to_use, 'all'))
%         species_to_use = 1:size(MRIO_threat_tensor, 1);
%     else
%         [~, ~, species_category_indexes] = intersect(species_characteristics.species_taxons, IUCN_threat_taxons, 'stable');
%         IUCN_threat_category = IUCN_threat_status(species_category_indexes);
%         species_to_use = find(ismember(IUCN_threat_category, analyse_MRIO_params.status_levels_to_use)); 
%     end
% 
%     if strcmp(analyse_MRIO_params.country_of_interest, 'global') 
%         industries_to_use = (1:length(industry_characteristics.country_names_list))';
%     else industries_to_use = find(strcmp(industry_characteristics.country_names_list, analyse_MRIO_params.country_of_interest));
%     end
%     
%     if strcmp(analyse_MRIO_params.industry_assessment_type, 'finalsale_based')
%         if strcmp(analyse_MRIO_params.analysis_type, 'by_country')
%             MRIO_subset = MRIO_threat_tensor(species_to_use, :, industries_to_use); 
%         else MRIO_subset = MRIO_threat_tensor(species_to_use, :, industries_to_use, :); 
%         end
%     else 
%         if strcmp(analyse_MRIO_params.analysis_type, 'by_country')
%             MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :);
%         else MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :, :);
%         end
%     
%     if strcmp(analyse_MRIO_params.country_of_interest, 'global') 
%         if strcmp(analyse_MRIO_params.analysis_type, 'by_country')
%             MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :);
%         else MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :, :);
%         end
%         
%         production_industries = industries_to_use(MRIO_subset.subs(:, 2));
%         final_sale_industries = industries_to_use(MRIO_subset.subs(:, 3));
%         
%     else
%         
%         if strcmp(analyse_MRIO_params.industry_assessment_type, 'production_based')
%             if strcmp(analyse_MRIO_params.analysis_type, 'by_country')
%                 MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :);
%             else MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :, :);
%             end
% 
%             production_industries = industries_to_use(MRIO_subset.subs(:, 2));
%             final_sale_industries = MRIO_subset.subs(:, 3);
%             
%         elseif strcmp(analyse_MRIO_params.industry_assessment_type, 'finalsale_based')
%             MRIO_subset = MRIO_threat_tensor(species_to_use, :, industries_to_use); 
%             production_industries = MRIO_subset.subs(:, 2);
%             final_sale_industries = industries_to_use(MRIO_subset.subs(:, 3));
%         end
%         
%     end
%     
%     MRIO_species_industry_identifiers = [final_sale_industries production_industries MRIO_subset.subs(:, 1)];
%     MRIO_species_threat_proportions = MRIO_subset.vals;
%     
%     if ~strcmp(analyse_MRIO_params.assessment_scale, 'global')
%         industries_to_use = industry_characteristics.country_index_list(final_sale_industries) - industry_characteristics.country_index_list(production_industries) == 0;
%         if strcmp(analyse_MRIO_params.assessment_scale, 'international')
%             industries_to_use = ~industries_to_use;
%         end
%         
%         MRIO_species_industry_identifiers = MRIO_species_industry_identifiers(industries_to_use, :);
%         MRIO_species_threat_proportions = MRIO_species_threat_proportions(industries_to_use);
%         
%     end
%         
%     
% end
% 


function [MRIO_species_industry_identifiers, MRIO_species_threat_proportions] = select_MRIO_subset(analyse_MRIO_params, MRIO_threat_tensor, species_characteristics, ...
                                                                                                industry_characteristics, IUCN_threat_taxons, IUCN_threat_status)
    
    if (strcmp(analyse_MRIO_params.status_levels_to_use, 'all'))
        species_to_use = 1:size(MRIO_threat_tensor, 1);
    else
        [~, ~, species_category_indexes] = intersect(species_characteristics.species_taxons, IUCN_threat_taxons, 'stable');
        IUCN_threat_category = IUCN_threat_status(species_category_indexes);
        species_to_use = find(ismember(IUCN_threat_category, analyse_MRIO_params.status_levels_to_use)); 
    end

    if strcmp(analyse_MRIO_params.country_of_interest, 'global') 
        industries_to_use = (1:length(industry_characteristics.country_names_list))';
        
        if numel(species_to_use) > 1
            MRIO_subset_spc = MRIO_threat_tensor(species_to_use, :, :);
        else
            MRIO_subset_spc = sptensor([],[], [1, size(MRIO_threat_tensor, 2) size(MRIO_threat_tensor, 3)]);
            set_to_use = MRIO_threat_tensor.subs(:, 1) == species_to_use;
            subs_to_use = MRIO_threat_tensor.subs(set_to_use, :, :);
            MRIO_subset_spc([ones(size(set_to_use)) subs_to_use(:, 2:3)]) = MRIO_threat_tensor.vals(set_to_use);
        end
        
        if numel(industries_to_use) > 1
            MRIO_subset = MRIO_subset_spc(:, industries_to_use, :);
        else
            MRIO_subset = sptensor([],[], [1, size(MRIO_subset_spc, 2) size(MRIO_subset_spc, 3)]);
            set_to_use = MRIO_subset_spc.subs(:, 2) == industries_to_use;
            subs_to_use = MRIO_subset_spc.subs(:, set_to_use, :);
            MRIO_subset([subs_to_use(:, 1) ones(size(set_to_use)) subs_to_use(:, 3)]) = MRIO_subset_spc.vals(set_to_use);
        end
        
        production_industries = industries_to_use(MRIO_subset.subs(:, 2));
        consumption_industries = industries_to_use(MRIO_subset.subs(:, 3));
        
    else
        disp('fix this routine set to above')
        
       
        industries_to_use = find(strcmp(industry_characteristics.country_names_list, analyse_MRIO_params.country_of_interest));
        if strcmp(analyse_MRIO_params.industry_assessment_type, 'production_based')
            MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :);
            
            production_industries = industries_to_use(MRIO_subset.subs(:, 2));
            consumption_industries = MRIO_subset.subs(:, 3);
            
        elseif strcmp(analyse_MRIO_params.industry_assessment_type, 'finalsale_based')
            MRIO_subset = MRIO_threat_tensor(species_to_use, :, industries_to_use); 
            production_industries = MRIO_subset.subs(:, 2);
            consumption_industries = industries_to_use(MRIO_subset.subs(:, 3));
        end
        
    end
    
    MRIO_species_industry_identifiers = [consumption_industries production_industries MRIO_subset.subs(:, 1)];
    MRIO_species_threat_proportions = MRIO_subset.vals;
    
    if ~strcmp(analyse_MRIO_params.assessment_scale, 'global')
        industries_to_use = industry_characteristics.country_index_list(consumption_industries) - industry_characteristics.country_index_list(production_industries) == 0;
        
        if strcmp(analyse_MRIO_params.assessment_scale, 'international')
            industries_to_use = ~industries_to_use;
        end
        
        MRIO_species_industry_identifiers = MRIO_species_industry_identifiers(industries_to_use, :);
        MRIO_species_threat_proportions = MRIO_species_threat_proportions(industries_to_use);
        
    end
        
    
end


function [species_threat_num] = build_species_threat_num(IUCN_species_names, unique_IUCN_species)
    species_threat_num = zeros(numel(unique_IUCN_species), 1);
    
    for i = 1:numel(unique_IUCN_species)
        species_threat_num(i) = sum(strcmp(IUCN_species_names, unique_IUCN_species(i)));
    end

end


function [consumption_industry_indexes, production_industry_indexes] = select_production_consumption_industry_indexes(sorted_MRIO_identifiers, analyse_MRIO_params, industries_to_use)
    if strcmp(analyse_MRIO_params.industry_assessment_type, 'production_based')
        production_industry_indexes = industries_to_use(sorted_MRIO_identifiers(:, analyse_MRIO_params.production_col));  % subset of industries means indexing goes from 1:N rather than the actual industry indexing
        consumption_industry_indexes = sorted_MRIO_identifiers(:, analyse_MRIO_params.consumption_col);
    
    elseif strcmp(analyse_MRIO_params.industry_assessment_type, 'finalsale_based')
        production_industry_indexes = sorted_MRIO_identifiers(:, analyse_MRIO_params.production_col);  
        consumption_industry_indexes = industries_to_use(sorted_MRIO_identifiers(:, analyse_MRIO_params.consumption_col));    
    end 
    
end

function sorted_MRIO_data = build_consumption_production_array(industry_characteristics, consumption_industry_indexes, production_industry_indexes)

    sorted_MRIO_data = cell(1, 4);
    sorted_MRIO_data{1} = industry_characteristics.country_names_list( consumption_industry_indexes ); 
    sorted_MRIO_data{2} = industry_characteristics.commodity_classification_list(consumption_industry_indexes);
    sorted_MRIO_data{3} = industry_characteristics.country_names_list( production_industry_indexes );  
    sorted_MRIO_data{4} = industry_characteristics.commodity_classification_list(production_industry_indexes);
    
end


% function assess_species_characteristics(species_characteristics, IUCN_data_object)
%     [~, ~, species_category_indexes] = intersect(species_characteristics.species_taxons, IUCN_data_object.IUCN_threat_taxons, 'stable');
%     T = cell2table([IUCN_data_object.IUCN_country_names num2cell(cellfun(@(x) sum(x), IUCN_data_object.x))], 'VariableNames', {'Country', 'GDP'})
%     
% end
% 
% function write_net_outputs(ranked_consumption_country_data, industry_characteristics, ranked_consumption_countries, ranked_production_country_data, ranked_production_countries, analyse_MRIO_params)
%         ranked_consumption_country_data_to_use = append_zero_block(ranked_consumption_country_data, industry_characteristics, ranked_consumption_countries);
%         ranked_production_country_data_to_use = append_zero_block(ranked_production_country, industry_characteristics, ranked_production_countries);
%         
%         [~, ~, sorted_production_inds] = intersect(ranked_consumption_country_data_to_use(:, 1), ranked_production_country_data_to_use(:, 1), 'stable');
%         
%         ranked_production_country_data_to_use = ranked_production_country_data_to_use(sorted_production_inds, :);
%         net_country_data = cell2mat(ranked_consumption_country_data(:, 2:end)) - cell2mat(ranked_production_country_data_to_use(:, 2:end));
%         
%         if (analyse_MRIO_params.sort_data == TRUE)
%         [~, sorted_country_indexes] = sort(net_country_data(:, end), 'ascend');
%                 net_country_data = [ranked_consumption_country_data(sorted_country_indexes, 1), num2cell(net_country_data(sorted_country_indexes, :))];
%         else net_country_data = [ranked_consumption_country_data(:, 1), num2cell(net_country_data)];
%         end
%             
% 
%         table_names = {'Country', 'Threatened_Animalia',  'Threatened_Plantae', 'Unique_Threatened_Species', 'Net_Threats'};
%         
%         filename = strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_net_', analyse_MRIO_params.assessment_scale, '_assessment_scale.txt');
%         T = cell2table(net_country_data, 'VariableNames', table_names);
%      
%         writetable(T, filename, 'Delimiter', 'tab')
% end


