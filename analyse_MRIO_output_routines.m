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

    [ranked_industry_paths_per_species, ranked_aggregated_threats_per_species, grouped_MRIO_species_aggregates, grouped_MRIO_industry_threat_vals] = aggregate_and_sort_paths(MRIO_species_threat_proportions, MRIO_species_industry_triplets(:, 1:2), MRIO_species_industry_triplets(:, 3));
   
    write_species_level_ranked_outputs(industry_characteristics, species_characteristics, analyse_MRIO_params, ranked_industry_paths_per_species, ranked_aggregated_threats_per_species, grouped_MRIO_species_aggregates)
    
    if analyse_MRIO_params.write_consumption_country_ranks == true
        [ranked_consumption_countries, ranked_consumption_threat_proportions_per_country, ranked_species_threatened_by_consumption_per_country, ranked_consumption_country_industry_threat_vals] = aggregate_and_sort_paths(ranked_aggregated_threats_per_species, industry_characteristics.country_index_list(ranked_industry_paths_per_species(:, 1)), grouped_MRIO_species_aggregates);
        ranked_consumption_country_data = build_ranked_outputs('by_country', 'consumption', analyse_MRIO_params, industry_characteristics, species_characteristics, ranked_consumption_countries, ... 
                                                                    ranked_consumption_threat_proportions_per_country, ranked_species_threatened_by_consumption_per_country);
    end
    
    if analyse_MRIO_params.write_production_country_ranks == true
        [ranked_production_countries, ranked_production_threat_proportions_per_country, ranked_species_threatened_by_production_per_country, ranked_production_country_industry_threat_vals] = aggregate_and_sort_paths(ranked_aggregated_threats_per_species, industry_characteristics.country_index_list(ranked_industry_paths_per_species(:, 2)), grouped_MRIO_species_aggregates);
        ranked_production_country_data = build_ranked_outputs('by_country', 'production', analyse_MRIO_params, industry_characteristics, species_characteristics, ranked_production_countries, ... 
                                                                    ranked_production_threat_proportions_per_country, ranked_species_threatened_by_production_per_country);
    end
    
    if analyse_MRIO_params.write_net_country_ranks == true
        write_net_outputs(ranked_consumption_country_data, industry_characteristics, ranked_consumption_countries, ranked_production_country_data, ranked_production_countries, analyse_MRIO_params)
    end
    
    if analyse_MRIO_params.write_industry_ranks == true
        [ranked_aggregated_industries_per_path, ranked_threat_proportions_per_industry, grouped_MRIO_industry_paths, grouped_aggregated_MRIO_industry_vals] = aggregate_and_sort_paths(ranked_aggregated_threats_per_species, ranked_industry_paths_per_species(:, 1), ranked_industry_paths_per_species(:, 2));

        [ranked_net_aggregated_industries, ranked_aggregated_threat_proportions_per_industry, industry_grouped_MRIO_species_aggregates, ~] = aggregate_and_sort_paths(ranked_aggregated_threats_per_species, ranked_industry_paths_per_species(:, 1), grouped_MRIO_species_aggregates);
   
        ranked_industry_data = build_ranked_outputs('by_industry', '', analyse_MRIO_params, industry_characteristics, species_characteristics, ranked_aggregated_industries_per_path, ranked_aggregated_threat_proportions_per_industry, industry_grouped_MRIO_species_aggregates);
               
        if analyse_MRIO_params.write_expanded_table == true
            expanded_MRIO_data = write_expanded_ranked_outputs(ranked_industry_data, analyse_MRIO_params, industry_characteristics, species_characteristics, grouped_MRIO_industry_paths, industry_grouped_MRIO_species_aggregates, grouped_aggregated_MRIO_industry_vals);
        end
        
    end
    
    
end

function assess_species_characteristics(species_characteristics, IUCN_data_object)
    [~, ~, species_category_indexes] = intersect(species_characteristics.species_taxons, IUCN_data_object.IUCN_threat_taxons, 'stable');
    T = cell2table([IUCN_data_object.IUCN_country_names num2cell(cellfun(@(x) sum(x), IUCN_data_object.x))], 'VariableNames', {'Country', 'GDP'})
    
end

function write_net_outputs(ranked_consumption_country_data, industry_characteristics, ranked_consumption_countries, ranked_production_country_data, ranked_production_countries, analyse_MRIO_params)
        ranked_consumption_country_data_to_use = append_zero_block(ranked_consumption_country_data, industry_characteristics, ranked_consumption_countries);
        ranked_production_country_data_to_use = append_zero_block(ranked_production_country, industry_characteristics, ranked_production_countries);
        
        [~, ~, sorted_production_inds] = intersect(ranked_consumption_country_data_to_use(:, 1), ranked_production_country_data_to_use(:, 1), 'stable');
        
        ranked_production_country_data_to_use = ranked_production_country_data_to_use(sorted_production_inds, :);
        net_country_data = cell2mat(ranked_consumption_country_data(:, 2:end)) - cell2mat(ranked_production_country_data_to_use(:, 2:end));
        
        [~, sorted_country_indexes] = sort(net_country_data(:, end), 'ascend');
        
        net_country_data = [ranked_consumption_country_data(sorted_country_indexes, 1), num2cell(net_country_data(sorted_country_indexes, :))];
        table_names = {'Country', 'Threatened_Animalia',  'Threatened_Plantae', 'Unique_Threatened_Species', 'Net_Threats'};
        
        filename = strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_net_', analyse_MRIO_params.assessment_scale, '_assessment_scale.txt');
        T = cell2table(net_country_data, 'VariableNames', table_names);
     
        writetable(T, filename, 'Delimiter', 'tab')
end

function write_species_level_ranked_outputs(industry_characteristics, species_characteristics, analyse_MRIO_params, ranked_industry_paths_per_species, ranked_aggregated_threats_per_species, grouped_MRIO_species_aggregates)
       
    ranked_MRIO_industries = cell(length(ranked_industry_paths_per_species), 4);
    ranked_MRIO_industries(:, 1) = industry_characteristics.country_names_list( ranked_industry_paths_per_species(:, 1)); 
    ranked_MRIO_industries(:, 2) = industry_characteristics.commodity_classification_list(ranked_industry_paths_per_species(:, 1));
    ranked_MRIO_industries(:, 3) = industry_characteristics.country_names_list( ranked_industry_paths_per_species(:, 2)); 
    ranked_MRIO_industries(:, 4) = industry_characteristics.commodity_classification_list(ranked_industry_paths_per_species(:, 2));
    
    decomposed_aggregate_species_characteristics = cell(length(grouped_MRIO_species_aggregates), 2);

    decomposed_aggregate_species_characteristics(:, 1) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'ANIMALIA')))), grouped_MRIO_species_aggregates));
    decomposed_aggregate_species_characteristics(:, 2) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'PLANTAE')))), grouped_MRIO_species_aggregates));
    
    %decomposed_aggregate_species_characteristics(:, 3) = num2cell(cellfun(@(x) length(unique(x)), grouped_MRIO_species_aggregates));
    %decomposed_aggregate_species_characteristics(:, 4) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'FUNGI')))), grouped_MRIO_species_aggregates));
    %decomposed_aggregate_species_characteristics(:, 5) = num2cell(cellfun(@(x) length(unique(x(strcmp(species_characteristics.species_kingdom(x), 'CHROMISTA')))), grouped_MRIO_species_aggregates));
    
    
    ranked_MRIO_industries = [ranked_MRIO_industries decomposed_aggregate_species_characteristics num2cell(round(ranked_aggregated_threats_per_species, 1))];
    
    filename = strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_species_level_', analyse_MRIO_params.assessment_scale, '_assessment_scale.txt');
    table_names = {'Consumption_Country', 'Consumption_Industry', 'Production_Country', 'Production_Industry', 'Threatened_Animalia',  'Threatened_Plantae', 'Aggregated_Threats'};
       
    T = cell2table(ranked_MRIO_industries, 'VariableNames', table_names);
     
    writetable(T, filename, 'Delimiter', 'tab')
    
end



function ranked_countries = append_zero_block(ranked_countries, industry_characteristics, country_indexes_used)
    consumption_countries_to_append = setdiff(cell2mat(industry_characteristics.country_index_map(:, 2)), country_indexes_used);
        
    cellblock_to_append = cell(length(consumption_countries_to_append), size(ranked_countries, 2));
    cellblock_to_append(:, 1) = industry_characteristics.country_index_map(consumption_countries_to_append, 1);
        
    cellblock_to_append(:, 2:end) = num2cell(zeros(size(cellblock_to_append, 1), (size(cellblock_to_append, 2) - 1)));
    
    ranked_countries = vertcat(ranked_countries, cellblock_to_append);
end

function decomposed_aggregate_species_characteristics = decompose_aggregated_species(industry_grouped_species_aggregates, species_characteristics)
    
    decomposed_aggregate_species_characteristics = cell(length(industry_grouped_species_aggregates), 2);
    decomposed_aggregate_species_characteristics(:, 1) = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'ANIMALIA')), x, 'un', 0) ))), industry_grouped_species_aggregates));
    decomposed_aggregate_species_characteristics(:, 2) = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'PLANTAE')), x, 'un', 0) ))), industry_grouped_species_aggregates));
%    decomposed_aggregate_species_characteristics(:, 3) = num2cell(cellfun(@(x) length(unique(cell2mat(x))), industry_grouped_species_aggregates));
%    decomposed_aggregate_species_characteristics(:, 3) = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'FUNGI')), x, 'un', 0) ))), industry_grouped_species_aggregates));
%    decomposed_aggregate_species_characteristics(:, 4) = num2cell(cellfun(@(x) length(unique(cell2mat(cellfun(@(y) y(strcmp(species_characteristics.species_kingdom(y), 'CHROMISTA')), x, 'un', 0) ))), industry_grouped_species_aggregates));
    
end

function ranked_MRIO_industries = build_ranked_outputs(rank_type, file_identifier, analyse_MRIO_params, industry_characteristics, species_characteristics, object_to_aggregate_over, net_threat_aggregate, industry_grouped_species_aggregates)
    
    if strcmp(rank_type, 'by_industry')
        ranked_MRIO_industries = cell(length(object_to_aggregate_over), 2);
        ranked_MRIO_industries(:, 1) = industry_characteristics.country_names_list( object_to_aggregate_over);
        ranked_MRIO_industries(:, 2) = industry_characteristics.commodity_classification_list(object_to_aggregate_over);
        table_names = {'Consumption_Country', 'Consumption_Industry', 'Threatened_Animalia',  'Threatened_Plantae', 'Aggregated_Threats'};
        filename = strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_', file_identifier, '_', analyse_MRIO_params.assessment_scale, '_industry_assessment_scale.txt');
    else 
        ranked_MRIO_industries = cell(length(object_to_aggregate_over), 1);
        ranked_MRIO_industries(:, 1) = industry_characteristics.country_index_map( object_to_aggregate_over, 1);
        table_names = {'Consumption_Country', 'Threatened_Animalia',  'Threatened_Plantae', 'Aggregated_Threats'};
        filename = strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_', file_identifier, '_', analyse_MRIO_params.assessment_scale, '_assessment_scale.txt');
    end
    
    decomposed_aggregate_species_characteristics = decompose_aggregated_species(industry_grouped_species_aggregates, species_characteristics);
    ranked_MRIO_industries = [ranked_MRIO_industries decomposed_aggregate_species_characteristics  num2cell(round(net_threat_aggregate, 1))];
    
    ranked_MRIO_industries_to_output = append_zero_block(ranked_MRIO_industries, industry_characteristics, object_to_aggregate_over);
        
    T = cell2table(ranked_MRIO_industries_to_output, 'VariableNames', table_names);
    writetable(T, filename, 'Delimiter', 'tab')

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
        writetable(T, strcat(analyse_MRIO_params.output_folder, analyse_MRIO_params.country_of_interest, '_', analyse_MRIO_params.assessment_type, '_', analyse_MRIO_params.assessment_scale, '_assessment_scale_expanded.txt'), 'Delimiter', 'tab')
    
end




function [aggregated_paths, aggregated_path_vals, grouped_aggregates, grouped_path_vals] = aggregate_and_sort_paths(vals_to_aggregate, paths_to_aggregate, objects_to_aggregate)
    
    [aggregated_paths, ~, path_indexes] = unique(paths_to_aggregate, 'rows', 'stable');
    
    aggregated_path_vals = accumarray(path_indexes, vals_to_aggregate);
    
    grouped_path_indexes = accumarray(path_indexes, find(path_indexes), [], @(rows){rows});  
    [aggregated_path_vals, sorted_indexes] = sort(aggregated_path_vals, 'descend');
    aggregated_paths = aggregated_paths(sorted_indexes, :);
    grouped_path_indexes = grouped_path_indexes(sorted_indexes);

    sorted_grouped_path_indexes = cellfun(@(x) x(sort_indexes(vals_to_aggregate(x))), grouped_path_indexes, 'un', false);
    grouped_path_vals = cellfun(@(x) vals_to_aggregate(x), sorted_grouped_path_indexes, 'un', false);
    grouped_aggregates = cellfun(@(x) objects_to_aggregate(x, :), sorted_grouped_path_indexes, 'un', false);
    
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
     
    if strcmp(analyse_MRIO_params.assessment_type, 'production_based')
        domestic_threat_indexes = strcmp(MRIO_consumption_production{1}, analyse_MRIO_params.country_of_interest);
    elseif strcmp(analyse_MRIO_params.assessment_type, 'consumption_based')
        domestic_threat_indexes = strcmp(MRIO_consumption_production{3}, analyse_MRIO_params.country_of_interest);
    end

    domestic_threat_paths = global_threat_paths(domestic_threat_indexes, :);
    international_threat_paths = global_threat_paths(~domestic_threat_indexes, :);
    
end
    
function save_MRIO_threat_rankings(domestic_threat_paths, international_threat_paths, global_threat_paths)
    if strcmp(analyse_MRIO_params.assessment_type, 'production_based')
        
        file_prefix = strcat(analyse_MRIO_params.output_folder, 'threat_proportion_rankings_consumption_type_');
        file_suffix = strcat('_production_country_', analyse_MRIO_params.country_of_interest, '.mat');
        int_filename = strcat(file_prefix, 'international', file_suffix);
        dom_filename = strcat(file_prefix, 'domestic', file_suffix);
        global_filename = strcat(file_prefix, 'global', file_suffix);
    elseif strcmp(analyse_MRIO_params.assessment_type, 'consumption_based')
        
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

function MRIO_threat_tensor = build_threat_tensor(analyse_MRIO_params)

%     load([analyse_MRIO_params.datapath 'SpThrSubs_domestic_final.mat'])
%     load([analyse_MRIO_params.datapath 'SpThrVals_domestic_final.mat'])
    load([analyse_MRIO_params.threat_indexes_filename])
    load([analyse_MRIO_params.threat_vals_filename])

    MRIO_threat_tensor = sptensor(double(SpThrSubs), double(SpThrVals), double(max(SpThrSubs)));
    save([analyse_MRIO_params.datapath analyse_MRIO_params.MRIO_threat_tensor_filename], 'MRIO_threat_tensor')
    
end

function [MRIO_species_industry_triplets, MRIO_species_threat_proportions] = select_MRIO_subset(analyse_MRIO_params, MRIO_threat_tensor, species_characteristics, ...
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
        MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :);
        production_industries = industries_to_use(MRIO_subset.subs(:, 2));
        consumption_industries = industries_to_use(MRIO_subset.subs(:, 3));
    else
        industries_to_use = find(strcmp(industry_characteristics.country_names_list, analyse_MRIO_params.country_of_interest));
        if strcmp(analyse_MRIO_params.assessment_type, 'production_based')
            MRIO_subset = MRIO_threat_tensor(species_to_use, industries_to_use, :);
            production_industries = industries_to_use(MRIO_subset.subs(:, 2));
            consumption_industries = MRIO_subset.subs(:, 3);
            
        elseif strcmp(analyse_MRIO_params.assessment_type, 'consumption_based')
            MRIO_subset = MRIO_threat_tensor(species_to_use, :, industries_to_use); 
            production_industries = MRIO_subset.subs(:, 2);
            consumption_industries = industries_to_use(MRIO_subset.subs(:, 3));
        end
        
    end
    
    MRIO_species_industry_triplets = [consumption_industries production_industries MRIO_subset.subs(:, 1)];
    MRIO_species_threat_proportions = MRIO_subset.vals;
    
    if ~strcmp(analyse_MRIO_params.assessment_scale, 'global')
        industries_to_use = industry_characteristics.country_index_list(consumption_industries) - industry_characteristics.country_index_list(production_industries) == 0;
        if strcmp(analyse_MRIO_params.assessment_scale, 'international')
            industries_to_use = ~industries_to_use;
        end
        
        MRIO_species_industry_triplets = MRIO_species_industry_triplets(industries_to_use, :);
        MRIO_species_threat_proportions = MRIO_species_threat_proportions(industries_to_use);
        
    end
        
    
end


function [species_threat_num] = build_species_threat_num(IUCN_species_names, unique_IUCN_species)
    species_threat_num = zeros(numel(unique_IUCN_species), 1);
    
    for i = 1:numel(unique_IUCN_species)
        species_threat_num(i) = sum(strcmp(IUCN_species_names, unique_IUCN_species(i)));
    end

end


function [consumption_industry_indexes, production_industry_indexes] = select_production_consumption_industry_indexes(sorted_MRIO_triplets, analyse_MRIO_params, industries_to_use)
    if strcmp(analyse_MRIO_params.assessment_type, 'production_based')
        production_industry_indexes = industries_to_use(sorted_MRIO_triplets(:, analyse_MRIO_params.production_col));  % subset of industries means indexing goes from 1:N rather than the actual industry indexing
        consumption_industry_indexes = sorted_MRIO_triplets(:, analyse_MRIO_params.consumption_col);
    
    elseif strcmp(analyse_MRIO_params.assessment_type, 'consumption_based')
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




