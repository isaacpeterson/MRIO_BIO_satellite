function MRIO_outputs = analyse_MRIO_output_routines(analyse_MRIO_params, IUCN_data_object)
    tic
    disp('loading tensor objects')
    if (analyse_MRIO_params.build_threat_tensor == true)
        MRIO_threat_tensor = build_threat_tensor(analyse_MRIO_params);
    else 
        load([analyse_MRIO_params.datapath analyse_MRIO_params.MRIO_threat_tensor_filename])
    end

    load([analyse_MRIO_params.datapath analyse_MRIO_params.satellite_species_characteristics_filename]);
    [~, ~, inds_to_use] = intersect(species_characteristics.species_taxons, IUCN_data_object.IUCN_threat_taxons, 'stable');
    species_characteristics.species_kingdom = IUCN_data_object.IUCN_species_kingdom(inds_to_use);

    disp(['processing tensors...'])
    
    industry_characteristics = build_industry_characteristics(analyse_MRIO_params);

    [MRIO_species_industry_triplets, MRIO_species_threat_proportions] = select_MRIO_subset(analyse_MRIO_params, MRIO_threat_tensor, species_characteristics, ...
                                                                            industry_characteristics, IUCN_data_object.IUCN_threat_taxons, IUCN_data_object.IUCN_threat_status);
                                                                
    MRIO_consumption_production = build_consumption_production_array(industry_characteristics, MRIO_species_industry_triplets(:, 1), MRIO_species_industry_triplets(:, 2));


    
    if strcmp(analyse_MRIO_params.rank_type, 'by_industry')
        [aggregated_MRIO_industry_threat_paths, aggregated_MRIO_threats, grouped_MRIO_species_aggregates, grouped_MRIO_industry_threat_vals] = aggregate_paths(MRIO_species_threat_proportions, MRIO_species_industry_triplets(:, 1:2), MRIO_species_industry_triplets(:, 3));
   
        [aggregated_MRIO_industry_paths, aggregated_MRIO_vals, grouped_MRIO_industry_paths, grouped_MRIO_industry_vals] = aggregate_paths(aggregated_MRIO_threats, aggregated_MRIO_industry_threat_paths(:, 1), aggregated_MRIO_industry_threat_paths(:, 2));
    
        [~, net_MRIO_threat_proportion, industry_grouped_MRIO_species_aggregates, ~] = aggregate_paths(aggregated_MRIO_threats, aggregated_MRIO_industry_threat_paths(:, 1), grouped_MRIO_species_aggregates);
    
        sort_MRIO_outputs_by_aggregated_industry_threats(analyse_MRIO_params, industry_characteristics, species_characteristics, aggregated_MRIO_industry_paths, grouped_MRIO_industry_paths, grouped_MRIO_species_aggregates, net_MRIO_threat_proportion, industry_grouped_MRIO_species_aggregates, grouped_MRIO_industry_vals);
    else 
            [domestic_threat_array, international_threat_array, global_threat_array] = build_MRIO_threat_rankings(MRIO_consumption_production, species_characteristics, analyse_MRIO_params, MRIO_species_threat_proportions, ...
                                                                                 MRIO_species_industry_triplets, IUCN_data_object.IUCN_threat_taxons);
    end
    toc
%     aggregated_MRIO_threats = aggregate_MRIO_industry_paths(aggregated_MRIO_threats);
%     
%     aggregated_MRIO_data = build_consumption_production_array(industry_characteristics, aggregated_MRIO_threats(:, 1), aggregated_MRIO_threats(:, 2));
%     aggregated_MRIO_data{5} = aggregated_MRIO_threats(:, 3);
%     aggregated_MRIO_data = [aggregated_MRIO_data{:}];
% 
%     save(strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.threat_direction, '_', analyse_MRIO_params.country_of_interest, '_ranked_by_aggregated_industry_paths.mat'), 'aggregated_MRIO_data', '-v7.3')
%     
%     MRIO_outputs = struct();
%     MRIO_outputs.aggregated_MRIO_data = aggregated_MRIO_data;
%     MRIO_outputs.domestic_threat_array = domestic_threat_array;
%     MRIO_outputs.international_threat_array = international_threat_array;
%     MRIO_outputs.global_threat_array = global_threat_array;
    
end


function sorted_MRIO_data = sort_MRIO_outputs_by_aggregated_country_threats(analyse_MRIO_params, industry_characteristics, species_characteristics, aggregated_MRIO_industry_paths, grouped_MRIO_industry_paths, grouped_MRIO_species_aggregates, ...
                                             net_MRIO_threat_proportion, industry_grouped_MRIO_species_aggregates, grouped_MRIO_industry_vals)
    
    
    total_MRIO_species_industry_paths = num2cell(cellfun(@(x) length(unique(cell2mat(x))), industry_grouped_MRIO_species_aggregates));
    
    total_unique_animalia = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'ANIMALIA')), x, 'un', 0) ))), industry_grouped_MRIO_species_aggregates));
    total_unique_plantae = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'PLANTAE')), x, 'un', 0) ))), industry_grouped_MRIO_species_aggregates));
    total_unique_fungi = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'FUNGI')), x, 'un', 0) ))), industry_grouped_MRIO_species_aggregates));
    total_unique_chromista = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'CHROMISTA')), x, 'un', 0) ))), industry_grouped_MRIO_species_aggregates));
    
    sorted_MRIO_data = cell(1, 16);
    sorted_MRIO_data{1} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), industry_characteristics.country_names_list( aggregated_MRIO_industry_paths), grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{2} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), industry_characteristics.commodity_classification_list(aggregated_MRIO_industry_paths), grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{3} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), num2cell(round(net_MRIO_threat_proportion, 1)), grouped_MRIO_industry_paths, 'un', 0); 
    sorted_MRIO_data{4} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_MRIO_species_industry_paths, grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{5} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_unique_animalia, grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{6} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_unique_plantae, grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{7} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_unique_fungi, grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{8} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_unique_chromista, grouped_MRIO_industry_paths, 'un', 0);
    
    summed_MRIO_species_industry_paths = cellfun(@(x) cellfun(@(y) length(y), x,'un',0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data{9} = cellfun(@(x) industry_characteristics.country_names_list(x), grouped_MRIO_industry_paths, 'un', false); 
    sorted_MRIO_data{10} = cellfun(@(x) industry_characteristics.commodity_classification_list(x), grouped_MRIO_industry_paths, 'un', false); 
    sorted_MRIO_data{11} = cellfun(@(x) num2cell(round(x, 1)), grouped_MRIO_industry_vals, 'un', 0); 
    sorted_MRIO_data{12} = summed_MRIO_species_industry_paths;
    
    sorted_MRIO_data{13} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'ANIMALIA'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data{14} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'PLANTAE'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data{15} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'FUNGI'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data{16} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'CHROMISTA'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data = cellfun(@(x) vertcat(x{:}), sorted_MRIO_data, 'un', 0);
    
    T = cell2table([sorted_MRIO_data{:}], 'VariableNames', {'Consumption_Country', 'Consumption_Industry', 'Total_Summed_Threats', 'Total_Unique_Threatened_Species', 'Total_unique_Animalia', ...
                                                            'Total_unique_Plantae', 'Total_unique_Fungi', 'Total_unique_Chromista', 'Production_Country', 'Production_Industry', ...
                                                            'Summed_Threats_Per_Tradepath', 'Threatened_Species_Per_Tradepath', 'Threatened_Animalia', 'Threatened_Plantae', 'Threatened_Fungi', 'Threatened_Chromista'});
    writetable(T, strcat('~/Documents/', analyse_MRIO_params.country_of_interest, '_', analyse_MRIO_params.threat_direction, '_MRIO.txt'), 'Delimiter', 'tab')
    
end



function sorted_MRIO_data = sort_MRIO_outputs_by_aggregated_industry_threats(analyse_MRIO_params, industry_characteristics, species_characteristics, aggregated_MRIO_industry_paths, grouped_MRIO_industry_paths, grouped_MRIO_species_aggregates, ...
                                             net_MRIO_threat_proportion, industry_grouped_MRIO_species_aggregates, grouped_MRIO_industry_vals)
    
    
    total_MRIO_species_industry_paths = num2cell(cellfun(@(x) length(unique(cell2mat(x))), industry_grouped_MRIO_species_aggregates));
    
    total_unique_animalia = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'ANIMALIA')), x, 'un', 0) ))), industry_grouped_MRIO_species_aggregates));
    total_unique_plantae = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'PLANTAE')), x, 'un', 0) ))), industry_grouped_MRIO_species_aggregates));
    total_unique_fungi = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'FUNGI')), x, 'un', 0) ))), industry_grouped_MRIO_species_aggregates));
    total_unique_chromista = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'CHROMISTA')), x, 'un', 0) ))), industry_grouped_MRIO_species_aggregates));
    
    sorted_MRIO_data = cell(1, 16);
    sorted_MRIO_data{1} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), industry_characteristics.country_names_list( aggregated_MRIO_industry_paths), grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{2} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), industry_characteristics.commodity_classification_list(aggregated_MRIO_industry_paths), grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{3} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), num2cell(round(net_MRIO_threat_proportion, 1)), grouped_MRIO_industry_paths, 'un', 0); 
    sorted_MRIO_data{4} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_MRIO_species_industry_paths, grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{5} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_unique_animalia, grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{6} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_unique_plantae, grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{7} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_unique_fungi, grouped_MRIO_industry_paths, 'un', 0);
    sorted_MRIO_data{8} = cellfun(@(x, y) vertcat(x, repmat({''}, (length(y) - 1), 1)), total_unique_chromista, grouped_MRIO_industry_paths, 'un', 0);
    
    summed_MRIO_species_industry_paths = cellfun(@(x) cellfun(@(y) length(y), x,'un',0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data{9} = cellfun(@(x) industry_characteristics.country_names_list(x), grouped_MRIO_industry_paths, 'un', false); 
    sorted_MRIO_data{10} = cellfun(@(x) industry_characteristics.commodity_classification_list(x), grouped_MRIO_industry_paths, 'un', false); 
    sorted_MRIO_data{11} = cellfun(@(x) num2cell(round(x, 1)), grouped_MRIO_industry_vals, 'un', 0); 
    sorted_MRIO_data{12} = summed_MRIO_species_industry_paths;
    
    sorted_MRIO_data{13} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'ANIMALIA'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data{14} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'PLANTAE'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data{15} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'FUNGI'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data{16} = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(y), 'CHROMISTA'))), x,'un', 0), industry_grouped_MRIO_species_aggregates,'un', 0);
    sorted_MRIO_data = cellfun(@(x) vertcat(x{:}), sorted_MRIO_data, 'un', 0);
    
    T = cell2table([sorted_MRIO_data{:}], 'VariableNames', {'Consumption_Country', 'Consumption_Industry', 'Total_Summed_Threats', 'Total_Unique_Threatened_Species', 'Total_unique_Animalia', ...
                                                            'Total_unique_Plantae', 'Total_unique_Fungi', 'Total_unique_Chromista', 'Production_Country', 'Production_Industry', ...
                                                            'Summed_Threats_Per_Tradepath', 'Threatened_Species_Per_Tradepath', 'Threatened_Animalia', 'Threatened_Plantae', 'Threatened_Fungi', 'Threatened_Chromista'});
    writetable(T, strcat('~/Documents/', analyse_MRIO_params.country_of_interest, '_', analyse_MRIO_params.threat_direction, '_MRIO.txt'), 'Delimiter', 'tab')
end

function sorted_MRIO_data = build_industry_path_array(industry_characteristics, consumption_industry_indexes, production_industry_indexes)
    
    [aggregated_MRIO_industry_threat_paths, aggregated_MRIO_threats, grouped_MRIO_species_aggregates, grouped_MRIO_industry_threat_vals]
    sorted_MRIO_data = cell(1, 6);
    sorted_MRIO_data{1} = industry_characteristics.country_names_list( aggregated_MRIO_industry_threat_paths(:, 1)); 
    sorted_MRIO_data{2} = industry_characteristics.commodity_classification_list(aggregated_MRIO_industry_threat_paths(:, 1));
    sorted_MRIO_data{3} = industry_characteristics.country_names_list( aggregated_MRIO_industry_threat_paths(:, 2)); 
    sorted_MRIO_data{4} = industry_characteristics.commodity_classification_list(aggregated_MRIO_industry_threat_paths(:, 2));
    sorted_MRIO_data{5} = num2cell(aggregated_MRIO_threats);
    sorted_MRIO_data{6} = num2cell(cellfun('length', grouped_MRIO_species_aggregates));
    a = cellfun(@(x) species_characteristics.species_names(x), grouped_MRIO_species_aggregates, 'un', false);
    
end

function [aggregated_paths, aggregated_path_vals, grouped_aggregates, grouped_path_vals] = aggregate_paths(vals_to_aggregate, paths_to_aggregate, species_to_aggregate)
    
    [aggregated_paths, ~, path_indexes] = unique(paths_to_aggregate, 'rows', 'stable');
    aggregated_path_vals = accumarray(path_indexes, vals_to_aggregate);
    grouped_path_indexes = accumarray(path_indexes, find(path_indexes), [], @(rows){rows});  
    [aggregated_path_vals, sorted_indexes] = sort(aggregated_path_vals, 'descend');
    aggregated_paths = aggregated_paths(sorted_indexes, :);
    grouped_path_indexes = grouped_path_indexes(sorted_indexes);

    sorted_grouped_path_indexes = cellfun(@(x) x(sort_indexes(vals_to_aggregate(x))), grouped_path_indexes, 'un', false);
    grouped_path_vals = cellfun(@(x) vals_to_aggregate(x), sorted_grouped_path_indexes, 'un', false);
    grouped_aggregates = cellfun(@(x) species_to_aggregate(x, :), sorted_grouped_path_indexes, 'un', false);
    
end

function [sorted_inds] = sort_indexes(arr)
    [~, sorted_inds] = sort(arr, 'descend');
end


function [domestic_threat_paths, international_threat_paths, global_threat_paths] = build_MRIO_threat_rankings(MRIO_consumption_production, species_characteristics, analyse_MRIO_params, ...
                                                                                                            MRIO_species_threat_proportions, MRIO_species_industry_triplets, IUCN_species_taxons)
     
    MRIO_threat_set = build_MRIO_threat_set(species_characteristics.species_taxons(MRIO_species_industry_triplets(:, 3)), IUCN_species_taxons);
    
    MRIO_consumption_production{5} = species_characteristics.species_names(MRIO_species_industry_triplets(:, 3));
    MRIO_consumption_production{6} = num2cell(MRIO_threat_set);
    MRIO_consumption_production{7} = num2cell(MRIO_species_threat_proportions);
    global_threat_paths = [MRIO_consumption_production{:}];
    [~, sorted_indexes] = sort(MRIO_species_threat_proportions, 'descend');
    global_threat_paths = global_threat_paths(sorted_indexes, :);
     
    if strcmp(analyse_MRIO_params.threat_direction, 'production_based')
        domestic_threat_indexes = strcmp(MRIO_consumption_production{1}, analyse_MRIO_params.country_of_interest);
    elseif strcmp(analyse_MRIO_params.threat_direction, 'consumption_based')
        domestic_threat_indexes = strcmp(MRIO_consumption_production{3}, analyse_MRIO_params.country_of_interest);
    end

    domestic_threat_paths = global_threat_paths(domestic_threat_indexes, :);
    international_threat_paths = global_threat_paths(~domestic_threat_indexes, :);
    
end
    
function save_MRIO_threat_rankings(domestic_threat_paths, international_threat_paths, global_threat_paths)
    if strcmp(analyse_MRIO_params.threat_direction, 'production_based')
        
        file_prefix = strcat(analyse_MRIO_params.output_folder, 'threat_proportion_rankings_consumption_type_');
        file_suffix = strcat('_production_country_', analyse_MRIO_params.country_of_interest, '.mat');
        int_filename = strcat(file_prefix, 'international', file_suffix);
        dom_filename = strcat(file_prefix, 'domestic', file_suffix);
        global_filename = strcat(file_prefix, 'global', file_suffix);
    elseif strcmp(analyse_MRIO_params.threat_direction, 'consumption_based')
        
        file_prefix = strcat(analyse_MRIO_params.output_folder, 'threat_proportion_rankings_consumption_country_', analyse_MRIO_params.country_of_interest, '_production_type_');
        file_suffix =  strcat('.mat');
        int_filename = strcat(file_prefix, 'international', file_suffix);
        dom_filename = strcat(file_prefix, 'domestic' , file_suffix);
        global_filename = strcat(file_prefix, 'global', file_suffix);
    end
             
    save(int_filename, 'international_threat_paths', '-v7.3')
    save(dom_filename, 'domestic_threat_paths', '-v7.3')
    save(global_filename, 'global_threat_paths', '-v7.3')
end

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
    industry_characteristics = struct();
    industry_characteristics.country_names_list = x_data{1};
    industry_characteristics.country_codes_list = x_data{2};
    industry_characteristics.commodity_classification_list = x_data{4};
end

function MRIO_threat_tensor = build_threat_tensor(analyse_MRIO_params)

    load([analyse_MRIO_params.datapath 'SpThrSubs_domestic_final.mat'])
    load([analyse_MRIO_params.datapath 'SpThrVals_domestic_final.mat'])

    MRIO_threat_tensor = sptensor(double(SpThrSubs), double(SpThrVals), double(max(SpThrSubs)));
    save([analyse_MRIO_params.datapath analyse_MRIO_params.MRIO_threat_tensor_filename], 'MRIO_threat_tensor')
    
end

function [MRIO_species_industry_triplets, MRIO_species_threat_proportions] = select_MRIO_subset(analyse_MRIO_params, MRIO_threat_tensor, species_characteristics, industry_characteristics, IUCN_threat_taxons, IUCN_threat_status)
    
    if (strcmp(analyse_MRIO_params.status_levels_to_use, 'all'))
        species_to_use = 1:size(MRIO_threat_tensor, 1);
    else
        [~, ~, species_category_indexes] = intersect(species_characteristics.species_taxons, IUCN_threat_taxons, 'stable');
        IUCN_threat_category = IUCN_threat_status(species_category_indexes);
        species_to_use = find(ismember(IUCN_threat_category, analyse_MRIO_params.status_levels_to_use)); 
    end

    if strcmp(analyse_MRIO_params.country_of_interest, 'all') 
        industries_to_use = (1:length(industry_characteristics.country_names_list))';
        MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :);
        production_industries = industries_to_use(MRIO_subset.subs(:, 2));
        consumption_industries = industries_to_use(MRIO_subset.subs(:, 3));
    else
        industries_to_use = find(strcmp(industry_characteristics.country_names_list, analyse_MRIO_params.country_of_interest));
        if strcmp(analyse_MRIO_params.threat_direction, 'production_based')
            MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :);
            production_industries = industries_to_use(MRIO_subset.subs(:, 2));
            consumption_industries = MRIO_subset.subs(:, 3);
            
        elseif strcmp(analyse_MRIO_params.threat_direction, 'consumption_based')
            MRIO_subset = MRIO_threat_tensor(species_to_use, :, industries_to_use); 
            production_industries = MRIO_subset.subs(:, 2);
            consumption_industries = industries_to_use(MRIO_subset.subs(:, 3));
        end
    end
    MRIO_species_industry_triplets = [consumption_industries production_industries MRIO_subset.subs(:, 1)];
    MRIO_species_threat_proportions = MRIO_subset.vals;
end


function [species_threat_num] = build_species_threat_num(IUCN_species_names, unique_IUCN_species)
    species_threat_num = zeros(numel(unique_IUCN_species), 1);
    
    for i = 1:numel(unique_IUCN_species)
        species_threat_num(i) = sum(strcmp(IUCN_species_names, unique_IUCN_species(i)));
    end

end


function [consumption_industry_indexes, production_industry_indexes] = select_production_consumption_industry_indexes(sorted_MRIO_triplets, analyse_MRIO_params, industries_to_use)
    if strcmp(analyse_MRIO_params.threat_direction, 'production_based')
        production_industry_indexes = industries_to_use(sorted_MRIO_triplets(:, analyse_MRIO_params.production_col));  % subset of industries means indexing goes from 1:N rather than the actual industry indexing
        consumption_industry_indexes = sorted_MRIO_triplets(:, analyse_MRIO_params.consumption_col);
    
    elseif strcmp(analyse_MRIO_params.threat_direction, 'consumption_based')
        production_industry_indexes = sorted_MRIO_triplets(:, analyse_MRIO_params.production_col);  
        consumption_industry_indexes = industries_to_use(sorted_MRIO_triplets(:, analyse_MRIO_params.consumption_col));    
    end 
    
end

function sorted_MRIO_data = build_consumption_production_array(industry_characteristics, consumption_industry_indexes, production_industry_indexes)

    sorted_MRIO_data = cell(1, 4);
    sorted_MRIO_data{1} = industry_characteristics.country_names_list( consumption_industry_indexes ); 
    sorted_MRIO_data{2} = industry_characteristics.commodity_classification_list(consumption_industry_indexes);
    sorted_MRIO_data{3} = industry_characteristics.country_names_list( production_industry_indexes );  
    sorted_MRIO_data{4} = industry_characteristics.commodity_classification_list(production_industry_indexes);
    
end




