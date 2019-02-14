function trade_characteristics = analyse_global_consumption_routines(IUCN_data_object, analyse_MRIO_params, species_characteristics)

    trade_characteristics = struct();
    trade_characteristics.finalsale_data = build_finalsale_level_data(analyse_MRIO_params, IUCN_data_object, species_characteristics);
    
    [trade_characteristics.country_scale.unique_countries, country_indexes_to_use] = setdiff(unique(IUCN_data_object.industry_characteristics.country_names_list, 'stable'), analyse_MRIO_params.countries_to_exclude, 'stable');    
    trade_characteristics.country_scale.finalsale_data_grouped_by_consumption_country = build_consumption_level_data(analyse_MRIO_params, IUCN_data_object, species_characteristics, country_indexes_to_use');

    trade_characteristics.sector_to_sector_scale = expand_sector_to_sector_finasale_data(country_indexes_to_use, trade_characteristics.country_scale.finalsale_data_grouped_by_consumption_country, IUCN_data_object);

    trade_characteristics.aggregated_sector_scale = expand_aggregated_sector_scale_data(trade_characteristics.country_scale.finalsale_data_grouped_by_consumption_country, ...
                                                                                        trade_characteristics.country_scale.unique_countries, country_indexes_to_use);
                                                                                    
    trade_characteristics.aggregated_sector_scale = link_consumption_to_aggregated_finalsale(trade_characteristics.aggregated_sector_scale, trade_characteristics.finalsale_data, IUCN_data_object, analyse_MRIO_params);
   
    trade_characteristics.country_scale.species_counts = build_country_scale_species_counts(trade_characteristics.aggregated_sector_scale, ...
                                                                                            trade_characteristics.finalsale_data, ...
                                                                                            country_indexes_to_use, ...
                                                                                            analyse_MRIO_params, ...
                                                                                            species_characteristics);
    
    trade_characteristics.country_scale.consumption_characteristics = run_consumption_production_assessment(trade_characteristics.sector_to_sector_scale.consumption_country_list, trade_characteristics.sector_to_sector_scale, ...
                                                                                                            trade_characteristics.country_scale.unique_countries, country_indexes_to_use, trade_characteristics.country_scale.species_counts, analyse_MRIO_params);
    
    trade_characteristics.country_scale.production_characteristics = run_consumption_production_assessment(trade_characteristics.sector_to_sector_scale.production_country_list, trade_characteristics.sector_to_sector_scale, ...
                                                                                                        trade_characteristics.country_scale.unique_countries, country_indexes_to_use, trade_characteristics.country_scale.species_counts, analyse_MRIO_params);
                                                                                       
    trade_characteristics.country_scale.import_characteristics = build_import_trade_characteristics(trade_characteristics, analyse_MRIO_params, trade_characteristics.country_scale.unique_countries, country_indexes_to_use);
    
    trade_characteristics.country_scale.net_trade_characteristics = run_country_scale_net_trade_assessment(trade_characteristics);     
    
    trade_characteristics.global_scale = build_global_scale_characteristics(trade_characteristics);
end

function global_scale_characteristics = build_global_scale_characteristics(trade_characteristics)

    global_scale_characteristics = struct();
    global_scale_characteristics.global_total = sum(trade_characteristics.sector_to_sector_scale.finalsale_production_threat_intensities);
    global_scale_characteristics.international_trade.total = sum(trade_characteristics.sector_to_sector_scale.finalsale_production_threat_intensities(trade_characteristics.sector_to_sector_scale.international_indexes));
    
    low_income_import_sectors = vertcat(trade_characteristics.country_scale.import_characteristics.low_income.import_indexes{:});
    global_scale_characteristics.international_trade.low_income = sum(trade_characteristics.sector_to_sector_scale.finalsale_production_threat_intensities(low_income_import_sectors));
    
    high_income_import_sectors = vertcat(trade_characteristics.country_scale.import_characteristics.high_income.import_indexes{:});
    global_scale_characteristics.international_trade.high_income = sum(trade_characteristics.sector_to_sector_scale.finalsale_production_threat_intensities(high_income_import_sectors));
    
    global_scale_characteristics.international_trade_proportion = global_scale_characteristics.international_trade.total/global_scale_characteristics.global_total;
    global_scale_characteristics.low_income_international_trade_proportion = global_scale_characteristics.international_trade.low_income/global_scale_characteristics.international_trade.total ;
    import_characteristics.low_income_indexes
end

function country_scale_consumption = build_consumption_level_data(analyse_MRIO_params, IUCN_data_object, species_characteristics, country_indexes_to_use)
    
    consumption_subs = load([analyse_MRIO_params.datapath '/PostExclNov/SpThrSubs_domestic_final.mat'] );
    consumption_vals = load([analyse_MRIO_params.datapath '/PostExclNov/SpThrVals_domestic_final.mat']);
    consumption_countries = load([analyse_MRIO_params.datapath '/PostExclNov/SpThrCnts_domestic_final.mat']);

    consumption_vals_to_keep = consumption_vals.SpThrVals > analyse_MRIO_params.data_threshold;
    consumption_vals.SpThrVals = consumption_vals.SpThrVals(consumption_vals_to_keep);
    consumption_subs.SpThrSubs = consumption_subs.SpThrSubs(consumption_vals_to_keep, :);
    consumption_countries.SpThrCountries = consumption_countries.SpThrCountries(consumption_vals_to_keep);
    
    country_num = length(country_indexes_to_use);
    
    country_scale_consumption = cell(country_num, 1);
    
    for country_index = 1:country_num
    
        if strcmp(analyse_MRIO_params.analysis_type, 'by_country')    
            current_country_set = find(consumption_countries.SpThrCountries == country_indexes_to_use(country_index));
        else
            current_country_set = 1:length(consumption_countries.SpThrCountries);
        end

        if (length(current_country_set) > 1)
            MRIO_threat_tensor = sptensor(double(consumption_subs.SpThrSubs(current_country_set, :)), double(consumption_vals.SpThrVals(current_country_set)), double(max(consumption_subs.SpThrSubs(current_country_set, :))));
            country_scale_consumption{country_index} = aggregate_MRIO_at_finalsale_routines(analyse_MRIO_params, IUCN_data_object, MRIO_threat_tensor, species_characteristics);
        end
    
    end
    
end

function finalsale_data = build_finalsale_level_data(analyse_MRIO_params, IUCN_data_object, species_characteristics)

    finalsale_subs = load([analyse_MRIO_params.datapath '/PostExclMarch/SpThrSubs_domestic_final.mat'] );
    finalsale_vals = load([analyse_MRIO_params.datapath '/PostExclMarch/SpThrVals_domestic_final.mat'] );

    % threshold data
    finalsale_vals_to_keep = (finalsale_vals.SpThrVals > analyse_MRIO_params.data_threshold);
    finalsale_vals.SpThrVals = finalsale_vals.SpThrVals(finalsale_vals_to_keep);
    finalsale_subs.SpThrSubs = finalsale_subs.SpThrSubs(finalsale_vals_to_keep, :);

    finalsale_threat_tensor = sptensor(double(finalsale_subs.SpThrSubs), double(finalsale_vals.SpThrVals), double(max(finalsale_subs.SpThrSubs)));
    finalsale_data = aggregate_MRIO_at_finalsale_routines(analyse_MRIO_params, IUCN_data_object, finalsale_threat_tensor, species_characteristics);
    
    tmp_exclusions = strcmp(IUCN_data_object.industry_characteristics.commodity_classification_list(finalsale_data.aggregated_sector_scale.aggregated_paths), 'Education, Health and Other Services');
                     
    finalsale_data.aggregated_sector_scale = structfun(@(x) x(~tmp_exclusions), finalsale_data.aggregated_sector_scale, 'UniformOutput', false);
    
    finalsale_data.aggregated_sector_scale.species_counts = count_species_groups(finalsale_data.aggregated_sector_scale.grouped_aggregates, species_characteristics, analyse_MRIO_params.groups_to_count);
    
end

function aggregated_sector_scale = expand_aggregated_sector_scale_data(consumption_country_set, unique_countries, country_indexes_to_use)
    
    aggregated_sector_scale = struct();
    %aggregated_sector_scale.consumption_country_index_list = cellfun(@(x, y) repmat({y}, [size(x.aggregated_sector_scale.aggregated_paths, 1) 1]), consumption_country_set, num2cell((1:length(unique_countries))'), 'UniformOutput', false);
   
    aggregated_sector_scale.consumption_country_index_list = cellfun(@(x, y) repmat({y}, [size(x.aggregated_sector_scale.aggregated_paths, 1) 1]), consumption_country_set, num2cell(country_indexes_to_use), 'UniformOutput', false);
    aggregated_sector_scale.consumption_country_index_list = cell2mat(vertcat(aggregated_sector_scale.consumption_country_index_list{:}));

    aggregated_sector_scale.consumption_country_list = cellfun(@(x, y) repmat({y}, [size(x.aggregated_sector_scale.aggregated_paths, 1) 1]), consumption_country_set, unique_countries, 'UniformOutput', false);
    aggregated_sector_scale.consumption_country_list = vertcat(aggregated_sector_scale.consumption_country_list{:});

    aggregated_sector_scale.finalsale_sector_list = cellfun(@(x) num2cell(x.aggregated_sector_scale.aggregated_paths), consumption_country_set, 'UniformOutput', false);
    aggregated_sector_scale.finalsale_sector_list = cell2mat(vertcat(aggregated_sector_scale.finalsale_sector_list{:}));

    aggregated_sector_scale.threat_intensities = cellfun(@(x) num2cell(x.aggregated_sector_scale.aggregated_path_vals), consumption_country_set, 'UniformOutput', false);
    aggregated_sector_scale.threat_intensities = cell2mat(vertcat(aggregated_sector_scale.threat_intensities{:}));

end


function aggregated_sector_scale = link_consumption_to_aggregated_finalsale(aggregated_sector_scale, finalsale_data, IUCN_data_object, analyse_MRIO_params)
    
    finalsale_link_indexes = arrayfun(@(x) find(aggregated_sector_scale.finalsale_sector_list == x), finalsale_data.aggregated_sector_scale.aggregated_paths, 'UniformOutput', false);

    non_empties = cellfun('length', finalsale_link_indexes) > 0;
    
    finalsale_link_indexes = finalsale_link_indexes(non_empties);

    aggregated_sector_scale.consumption_countries_per_finalsale_commodity = cellfun(@(x) aggregated_sector_scale.consumption_country_list(x), ...
                                                                                finalsale_link_indexes, 'UniformOutput', false);

    aggregated_sector_scale.expanded_finalsale_commodity_names = cellfun(@(x) IUCN_data_object.industry_characteristics.commodity_classification_list(aggregated_sector_scale.finalsale_sector_list(x)), ...
                                                                        finalsale_link_indexes, 'UniformOutput', false);

    aggregated_sector_scale.expanded_finalsale_country_names = cellfun(@(x) IUCN_data_object.industry_characteristics.country_names_list(aggregated_sector_scale.finalsale_sector_list(x)), ...
                                                                        finalsale_link_indexes, 'UniformOutput', false);    
                                                                    
    aggregated_sector_scale.scaled_threat_proportion = cellfun(@(x, y) aggregated_sector_scale.threat_intensities(x) ./ (sum(aggregated_sector_scale.threat_intensities(x))), ...
                                                               finalsale_link_indexes, 'UniformOutput', false);

    aggregated_sector_scale.scaled_threat_intensity = cellfun(@(x, y) x.*y, aggregated_sector_scale.scaled_threat_proportion, num2cell(finalsale_data.aggregated_sector_scale.aggregated_path_vals(non_empties)), 'UniformOutput', false);
    
    aggregated_sector_scale.scaled_threat_intensity = vertcat(aggregated_sector_scale.scaled_threat_intensity{:});
    
    aggregated_sector_scale.expanded_net_threat_intensities = cellfun(@(x, y) num2cell(repmat(finalsale_data.aggregated_sector_scale.aggregated_path_vals(y), [size(x, 1) 1])), ...
                                                                finalsale_link_indexes, num2cell(1:length(finalsale_link_indexes))', 'UniformOutput', false);
    
    expanded_species_counts = cell2struct(cell(length(analyse_MRIO_params.groups_to_count), 1), analyse_MRIO_params.groups_to_count);                                       
    
    for group_ind = 1:length(analyse_MRIO_params.groups_to_count)
        expanded_species_counts.(analyse_MRIO_params.groups_to_count{group_ind}) = cellfun(@(x, y) num2cell(repmat(finalsale_data.aggregated_sector_scale.species_counts.(analyse_MRIO_params.groups_to_count{group_ind})(y), [size(x, 1) 1])), ...
                                                                                            finalsale_link_indexes, num2cell(1:length(finalsale_link_indexes))', 'UniformOutput', false); 
    end  
    
    aggregated_sector_scale.consumption_to_finalsale_block = [vertcat(aggregated_sector_scale.consumption_countries_per_finalsale_commodity{:})...
                                                              vertcat(aggregated_sector_scale.expanded_finalsale_country_names{:})...
                                                              vertcat(aggregated_sector_scale.expanded_finalsale_commodity_names{:})...
                                                              vertcat(expanded_species_counts.(analyse_MRIO_params.groups_to_count{1}){:}) ...
                                                              vertcat(expanded_species_counts.(analyse_MRIO_params.groups_to_count{2}){:}) ...
                                                              num2cell(aggregated_sector_scale.scaled_threat_intensity), ...
                                                              num2cell(vertcat(aggregated_sector_scale.scaled_threat_proportion{:})) ...
                                                              vertcat(aggregated_sector_scale.expanded_net_threat_intensities{:})];
 
    [~, indexes_sorted_by_threat_intensity] = sort(aggregated_sector_scale.scaled_threat_intensity, 'descend');
    
    aggregated_sector_scale.consumption_to_finalsale_table = cell2table([num2cell(1:length(aggregated_sector_scale.scaled_threat_intensity))'...
                                                                        aggregated_sector_scale.consumption_to_finalsale_block(indexes_sorted_by_threat_intensity, :)], ...
                                                             'VariableNames', [{'global_rank', 'consumption_country', 'finalsale_country', 'finalsale_industry'}, analyse_MRIO_params.groups_to_count, ...
                                                                               {'attributed_threat_intensity', 'finalsale_threat_intensity_proportion', 'net_threat_intensity'}]);
                                                              
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


function species_counts = build_country_scale_species_counts(aggregated_sector_scale, finalsale_data, country_indexes_to_use, analyse_MRIO_params, species_characteristics)

    industry_link_indexes = arrayfun(@(x) find(aggregated_sector_scale.finalsale_sector_list == x), finalsale_data.aggregated_sector_scale.aggregated_paths, 'UniformOutput', false);

    %industry_link_indexes = arrayfun(@(x) find(aggregated_sector_scale.finalsale_sector_list == x), finalsale_data.aggregated_sector_scale.aggregated_paths, 'UniformOutput', false);
    
    non_empties = cellfun('length', industry_link_indexes) > 0;
    industry_link_indexes_to_use = industry_link_indexes(non_empties);

    %%%%%%%%%%% run analysis on species
    finalsale_data_to_use = finalsale_data.aggregated_sector_scale.grouped_aggregates(non_empties); 

    finalsale_blocks = cellfun(@(x, y) repmat({vertcat(y{:})}, [size(x, 1) 1]), industry_link_indexes_to_use, finalsale_data_to_use, 'UniformOutput', false);
    finalsale_blocks = vertcat(finalsale_blocks{:});
     
    grouped_industry_index_links = cellfun(@(x) [aggregated_sector_scale.consumption_country_index_list(x) repmat(aggregated_sector_scale.finalsale_sector_list(x), [1, 2])], industry_link_indexes_to_use, 'UniformOutput', false);

    grouped_industry_index_links = vertcat(grouped_industry_index_links{:});

    country_blocks = arrayfun(@(x) finalsale_blocks(grouped_industry_index_links(:, 1) == x), country_indexes_to_use, 'un', false);
    country_blocks = cellfun(@(x) vertcat(x{:}), country_blocks, 'un', false);
    country_blocks = cellfun(@(x) unique(x), country_blocks, 'un', false);

    species_counts = cellfun(@(x) cellfun(@(y) length(find(strcmp(species_characteristics.species_kingdom(x), y))), analyse_MRIO_params.groups_to_count), country_blocks, 'un', false);
    species_counts = vertcat(species_counts{:});
    
end



function sector_to_sector_scale = expand_sector_to_sector_finasale_data(country_indexes_to_use, consumption_country_set, IUCN_data_object)
    
    sector_to_sector_scale = struct();
    sector_to_sector_scale.finalsale_production_sectors = cellfun(@(x) num2cell(x.sector_to_sector_scale.aggregated_paths), consumption_country_set, 'UniformOutput', false); 
    sector_to_sector_scale.finalsale_production_threat_intensities = cellfun(@(x) num2cell(x.sector_to_sector_scale.aggregated_path_vals), consumption_country_set, 'UniformOutput', false);
    sector_to_sector_scale.consumption_country_list = cellfun(@(x, y) repmat(y, [size(x, 1) 1]), ...
                                                              sector_to_sector_scale.finalsale_production_sectors, num2cell(country_indexes_to_use),'UniformOutput', false);
   
    sector_to_sector_scale.consumption_country_list = vertcat(sector_to_sector_scale.consumption_country_list{:});
    sector_to_sector_scale.finalsale_production_sectors = cell2mat(vertcat(sector_to_sector_scale.finalsale_production_sectors{:}));
    sector_to_sector_scale.finalsale_production_threat_intensities = cell2mat(vertcat(sector_to_sector_scale.finalsale_production_threat_intensities{:}));
    
    sector_to_sector_scale.finalsale_country_list = IUCN_data_object.industry_characteristics.country_index_list(sector_to_sector_scale.finalsale_production_sectors(:, 1));
    sector_to_sector_scale.production_country_list = IUCN_data_object.industry_characteristics.country_index_list(sector_to_sector_scale.finalsale_production_sectors(:, 2));
    sector_to_sector_scale.international_indexes_including_finalsale = sum( abs(diff([sector_to_sector_scale.consumption_country_list sector_to_sector_scale.finalsale_country_list sector_to_sector_scale.production_country_list], 1, 2)), 2) > 0;
    sector_to_sector_scale.international_indexes = sum( abs(diff([sector_to_sector_scale.consumption_country_list sector_to_sector_scale.production_country_list], 1, 2)), 2) > 0;
end


function country_scale_trade_characteristics = run_country_scale_net_trade_assessment(trade_characteristics)
 
    country_scale_trade_characteristics = struct();
    country_scale_trade_characteristics.import_export_discriminator = trade_characteristics.country_scale.consumption_characteristics.threat_intensities_international ...
                                                - trade_characteristics.country_scale.production_characteristics.threat_intensities_international;   
                                            
    country_scale_trade_characteristics.net_threat_intensity = trade_characteristics.country_scale.consumption_characteristics.threat_intensities_global ...
                                            - trade_characteristics.country_scale.production_characteristics.threat_intensities_international;
    
    country_scale_trade_characteristics.table = [trade_characteristics.country_scale.consumption_characteristics.table(:, 1:(end - 1))...
                                                 cell2table( num2cell([trade_characteristics.country_scale.import_characteristics.low_income.imported_threats...
                                                                       trade_characteristics.country_scale.import_characteristics.low_income.imported_proportion...
                                                                       trade_characteristics.country_scale.production_characteristics.threat_intensities_international...
                                                                       trade_characteristics.country_scale.production_characteristics.international_trade_proportion]),...
                                                         'VariableNames', {'low_income_imports', 'low_income_imported_proportion', 'exports', 'proportion_exported'})...
                                                 trade_characteristics.country_scale.consumption_characteristics.table(:, end), ...
                                                 cell2table(num2cell(country_scale_trade_characteristics.net_threat_intensity), 'VariableNames', {'net'})];
    
end



function country_scale_trade_characteristics = run_consumption_production_assessment(block_to_use, sector_to_sector_scale, unique_countries, country_indexes_to_use, species_counts, analyse_MRIO_params)
    
    country_scale_trade_characteristics = struct();
    country_scale_trade_characteristics.global_blocks = cellfun(@(x) block_to_use == x, num2cell(country_indexes_to_use), 'UniformOutput', false);
    country_scale_trade_characteristics.international_blocks = cellfun(@(x) x & sector_to_sector_scale.international_indexes, country_scale_trade_characteristics.global_blocks, 'UniformOutput', false);
    country_scale_trade_characteristics.threat_intensities_international = cellfun(@(x) sum(sector_to_sector_scale.finalsale_production_threat_intensities(x)), country_scale_trade_characteristics.international_blocks, 'UniformOutput', false);
    country_scale_trade_characteristics.threat_intensities_international = vertcat(country_scale_trade_characteristics.threat_intensities_international{:});
    
    country_scale_trade_characteristics.threat_intensities_global = cellfun(@(x) sum(sector_to_sector_scale.finalsale_production_threat_intensities(x)), country_scale_trade_characteristics.global_blocks, 'UniformOutput', false);
    country_scale_trade_characteristics.threat_intensities_global = vertcat(country_scale_trade_characteristics.threat_intensities_global{:});
    
    country_scale_trade_characteristics.threat_intensities_domestic = country_scale_trade_characteristics.threat_intensities_global - country_scale_trade_characteristics.threat_intensities_international;
    country_scale_trade_characteristics.international_trade_proportion = country_scale_trade_characteristics.threat_intensities_international./country_scale_trade_characteristics.threat_intensities_global;
        
    trade_characteristics_block = [species_counts ...
                                   country_scale_trade_characteristics.threat_intensities_domestic ...
                                   country_scale_trade_characteristics.threat_intensities_international ...
                                   country_scale_trade_characteristics.international_trade_proportion...
                                   country_scale_trade_characteristics.threat_intensities_global];
                                  
    country_scale_trade_characteristics.table = cell2table([unique_countries num2cell(trade_characteristics_block)], ...
        'VariableNames', [{'Country'}, analyse_MRIO_params.groups_to_count, {'domestic_trade', 'international_trade', 'international_proportion', 'global_trade'}]);
    
end


function import_characteristics = build_import_trade_characteristics(trade_characteristics, analyse_MRIO_params, unique_countries, country_indexes_to_use)
    
    import_characteristics = struct();
    
    fid = fopen(analyse_MRIO_params.low_income_countries_filename);
        import_characteristics.developing_countries = textscan(fid,'%s', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);

    [~, import_characteristics.low_income_indexes] = intersect(unique_countries, import_characteristics.developing_countries{1});
    
    international_trade_indexes = cellfun(@(x) find(x), trade_characteristics.country_scale.consumption_characteristics.international_blocks, 'UniformOutput', false);
    
    low_income_countries = country_indexes_to_use(import_characteristics.low_income_indexes);
    high_income_countries = setdiff(country_indexes_to_use, country_indexes_to_use(import_characteristics.low_income_indexes));
    
    import_characteristics.low_income = assess_imports(trade_characteristics, international_trade_indexes, country_indexes_to_use(low_income_countries), unique_countries);
    import_characteristics.high_income = assess_imports(trade_characteristics, international_trade_indexes, high_income_countries, unique_countries);
    
end

function import_characteristics = assess_imports(trade_characteristics, international_trade_indexes, import_country_indexes, unique_countries)
    
    import_characteristics = struct();
    
    import_characteristics.import_country_indexes = import_country_indexes;
    
    import_characteristics.import_indexes = cellfun(@(x) ismember(trade_characteristics.sector_to_sector_scale.production_country_list(x), import_characteristics.import_country_indexes), ...
                                                               international_trade_indexes, 'UniformOutput', false);
                                                           
    import_characteristics.import_indexes = cellfun(@(x, y) x(y), international_trade_indexes, import_characteristics.import_indexes, 'UniformOutput', false);
    
    import_characteristics.imported_threats = cellfun(@(x) sum(trade_characteristics.sector_to_sector_scale.finalsale_production_threat_intensities(x)), ...
                                                            import_characteristics.import_indexes);
                                  
    import_characteristics.imported_proportion = import_characteristics.imported_threats ./ trade_characteristics.country_scale.consumption_characteristics.threat_intensities_international;
    
    import_characteristics.table = cell2table([unique_countries num2cell([import_characteristics.imported_threats ...
                                                                                    import_characteristics.imported_proportion...
                                                                                    trade_characteristics.country_scale.consumption_characteristics.threat_intensities_global])],...
                                            'VariableNames', {'consumption_country', 'imports', 'proportional_imports', 'net_consumption'});
end

