function [satellite_object] = build_iucn_satellite(iucn_data_object, satellite_params)
    
    display(['building satellite to ' satellite_params.system_type ' specification'])
    
    satellite_object = struct();
    
    satellite_object.species_characteristics = build_species_characteristics(iucn_data_object, satellite_params);
    
    satellite_params = process_params(satellite_params, iucn_data_object, satellite_object.species_characteristics);
    
    if (satellite_params.build_domestic_satellite == true)    
        [satellite_object.domestic_satellite, satellite_object.satellite_params] = build_domestic_satellite(iucn_data_object, satellite_params, satellite_object.species_characteristics); 
    end
    
    if (satellite_params.build_global_satellite == true)
        satellite_object.global_satellite = build_global_satellite(iucn_data_object, satellite_params, satellite_object.species_characteristics);  
    end
    
    if ~strcmp(satellite_params.satellite_collapse_type, 'none')
        satellite_object.collapsed_satellite = run_collapse_satellite_routines(satellite_params, iucn_data_object, species_characteristics);
    end

    if satellite_params.display_satellite == true
        display_satellite(satellite_object, satellite_object.satellite_params)
    end

    display(['satellite build finished at ', num2str(toc)])
    
end

function satellite_object = run_collapse_satellite_routines(satellite_params, iucn_data_object, species_characteristics)
    
    satellite_object = struct();
    satellite_object.species_taxons = species_characteristics.sorted_species_taxons;
    satellite_object.sorted_species_classification = species_characteristics.sorted_species_classification;
    satellite_object.iucn_country_code_names = iucn_data_object.iucn_country_code_names;
    satellite_object.country_names = iucn_data_object.iucn_country_names;
    
    if strcmp(satellite_params.satellite_collapse_type, 'eora25')
        satellite_concordance = csvread([satellite_params.satellite_collapse_concordance_filename]);
        satellite_concordance = satellite_concordance';
        sector_classification_num = repmat(25, [iucn_data_object.ncoun, 1]);
    end
    
    if satellite_params.return_satellite == true
        satellite_object.collapsed_satellite = arrayfun(@(x) zeros(satellite_params.classification_num, x), sector_classification_num, 'UniformOutput', false);  
    end
    
    disp(['collapsing to ',satellite_params.satellite_collapse_type, '...'])
    for country_index = 1:iucn_data_object.ncoun
        
        if strcmp(satellite_params.output_file_type, 'mat')
        	current_country = load([satellite_params.satellite_filepath, iucn_data_object.iucn_country_code_names{country_index}, '_agg.mat']);
            fname = fieldnames(current_country);
            current_country = current_country.(fname{1});
            
        elseif strcmp(satellite_params.output_file_type, 'csv') 
            csvread([satellite_params.satellite_filepath, iucn_data_object.iucn_country_code_names{country_index}, '_satellite.csv'])
        end
               
        if satellite_params.return_satellite == true
        	satellite_object.collapsed_satellite{country_index} = current_country*satellite_concordance;
        end
        
        disp([iucn_data_object.iucn_country_code_names{country_index}, ' collapsed at ' num2str(toc)])
        
    end
    
end


function collapsed_tensor = safe_collapse(current_tensor, collapse_inds, dim_to_collapse)

    if numel(collapse_inds) == 1
        collapsed_tensor = current_tensor;
    else collapsed_tensor = collapse(current_tensor, dim_to_collapse);
    end
    
end


function collapsed_species_tensor = collapse_species(iucn_tensor, satellite_params)
    
    sz = size(iucn_tensor);
    unique_species_inds = satellite_params.unique_species_inds;
    collapsed_species_tensor = sptensor([sz(1) satellite_params.classification_num sz(3:end)]);
    
    for collapse_ind = 1:satellite_params.classification_num
        if collapse_ind < satellite_params.classification_num
            species_set_to_use = unique_species_inds(collapse_ind):unique_species_inds(collapse_ind + 1);
            current_collapsed_tensor = safe_collapse(iucn_tensor(:, species_set_to_use, :), species_set_to_use, 2);
        else 
            species_set_to_use = unique_species_inds(collapse_ind):size(iucn_tensor,1);
            current_collapsed_tensor = safe_collapse(iucn_tensor(:, species_set_to_use, :), species_set_to_use, 2);
        end
        collapsed_species_tensor(:, collapse_ind, :) = current_collapsed_tensor;
    end  

end


function [current_satellite, satellite_params] = sort_satellite(current_satellite, iucn_data_object, satellite_params)
                                                                                                                                                
    if strcmp(satellite_params.country_sort_type, 'none') 
        satellite_params.sorted_country_names = iucn_data_object.iucn_country_names;
    else
        if strcmp(satellite_params.country_sort_type, 'eora')
            [eora_industry_specification, eora_country_names, eora_mapping_vector] = find_eora_industry_specification(satellite_params.eora_x_filename, iucn_data_object.un_to_iucn_codes, iucn_data_object.iucn_country_code_names); 
            disp('sorting to eora specification');
            [current_satellite, satellite_params] = sort_satellite_to_eora(current_satellite, satellite_params, eora_industry_specification, eora_mapping_vector, eora_country_names, iucn_data_object.iucn_country_names,  satellite_params.classification_num);
        else
            [current_satellite, satellite_params] = sort_satellite_by_country(current_satellite, iucn_data_object, satellite_params);
        end
    end
    
end


function satellite_params = process_params(satellite_params, iucn_data_object, species_characteristics)

    if strcmp(satellite_params.domestic_threats_to_aggregate, 'all')
        satellite_params.domestic_threats_to_aggregate = find(~iucn_data_object.allcountriesflag)';
    end
    
    if strcmp(satellite_params.global_threats_to_aggregate, 'all')
        satellite_params.global_threats_to_aggregate = find(iucn_data_object.allcountriesflag)';
    end   
    
    if  strcmp(satellite_params.status_levels_to_use, 'all')
       satellite_params.status_inds_to_use = 1:length(iucn_data_object.iucn_status_names);
    else
        [~, satellite_params.status_inds_to_use] = intersect(iucn_data_object.iucn_status_names, satellite_params.status_levels_to_use);
    end

    if satellite_params.collapse_through_species_sort_type == true
        [~, unique_species_inds] = unique(species_characteristics.sorted_species_classification, 'stable');
        satellite_params.classification_num = numel(unique_species_inds);
    else satellite_params.classification_num = numel(species_characteristics.sorted_q_row);
    end
    
    satellite_params.sector_lengths = cellfun('length', iucn_data_object.x);
    
end

                                                         
function [eora_industry_specification, eora_country_names, eora_mapping_vector] = find_eora_industry_specification(eora_x_filename, un_to_iucn_codes, iucn_country_code_names)
    
    fid = fopen(eora_x_filename);
        x_data = textscan(fid,'%s %s %s %s %f', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);
 
    eora_country_names_list = x_data{1};
    eora_country_codes_list = x_data{2};
    industry_commodity_classification_list = x_data{3};
    
    [unique_eora_country_codes, unique_eora_inds] = unique(eora_country_codes_list, 'stable');

    [~, un_inds, iucn_inds] = intersect(un_to_iucn_codes.iucn_country_codes, iucn_country_code_names, 'stable');

    un_to_iucn_inds = zeros(length(un_to_iucn_codes.un_country_codes), 1);
    un_to_iucn_inds(un_inds) = iucn_inds;

    [~, ~, un_inds_to_use] = intersect(unique_eora_country_codes, un_to_iucn_codes.un_country_codes, 'stable');
    
    eora_mapping_vector = un_to_iucn_inds(un_inds_to_use);

    eora_country_names = eora_country_names_list(unique_eora_inds);

    eora_industry_specification = zeros(numel(eora_mapping_vector), 2);
    
    for country_ind = 1:numel(eora_mapping_vector)
        current_set = strcmp(eora_country_codes_list, unique_eora_country_codes(country_ind));   
        industries_list = find(current_set & strcmp(industry_commodity_classification_list, 'Industries'));
        commodities_list = find(current_set & strcmp(industry_commodity_classification_list, 'Commodities'));
        eora_industry_specification(country_ind, :) = [numel(industries_list), numel(commodities_list)];
    end
    
%     for dummy_ind = 1:length(mapping_vector)
%         current_set = find(strcmp(x_data{2}, dummy_countries(dummy_ind)));
%         industry_data_list = strcmp(x_data{3}(current_set), 'Industries');
%         
%         if isempty(find(industry_data_list))
%             industry_data_list = strcmp(x_data{3}(current_set), 'Commodities');
%         end
%         
%         dummy_satellite{dummy_ind} = zeros(classification_num, length(find(industry_data_list)));
% 
%     end
  
end
    
        
    
function dummy_satellite = build_dummy_satellite(mapping_vector, sector_matches, classification_num, eora_satellite_params)

    dummy_satellite = cell(1, length(mapping_vector));
    
    for country_ind = 1:length(mapping_vector)
        if (mapping_vector(country_ind) > 0) && (sector_matches(country_ind) == true)
            dummy_satellite{country_ind} = zeros(classification_num, eora_satellite_params(country_ind, 2));   
        elseif (mapping_vector(country_ind) == 0) 
             dummy_satellite{country_ind} = zeros(classification_num, sum(eora_satellite_params(country_ind, :)));
        end        
    end
    
end

function [expanded_satellite, satellite_params] = sort_satellite_to_eora(current_satellite, satellite_params, eora_industry_specification, eora_mapping_vector, eora_country_names, iucn_country_names, classification_num) 
    
    dummy_inds = find(eora_mapping_vector == 0); 
    dummy_country_names = eora_country_names(dummy_inds);    
    expanded_satellite = cell(1, length(eora_mapping_vector));
    expanded_satellite((eora_mapping_vector > 0)) = current_satellite(eora_mapping_vector(eora_mapping_vector > 0));
    sector_lengths = cellfun(@(x) size(x, 2), expanded_satellite);
    sector_matches = (sector_lengths' == eora_industry_specification(:, 1));
        
    dummy_satellite = build_dummy_satellite(eora_mapping_vector, sector_matches, classification_num, eora_industry_specification);
    
    expanded_satellite = cellfun(@(x, y) [x y], expanded_satellite, dummy_satellite, 'UniformOutput', false);
    
    satellite_params.sorted_country_names(dummy_inds) = dummy_country_names;
    satellite_params.sorted_country_names(eora_mapping_vector > 0) = iucn_country_names(eora_mapping_vector(eora_mapping_vector > 0));

    satellite_params.sector_lengths = cellfun(@(x) size(x, 2), expanded_satellite);
    
end
    

function check_satellite(split_satellite, iucn_data_object, status_inds_to_use, threat_type)
    industry_threat_counts = cellfun(@(x) sum(x, 2), split_satellite, 'UniformOutput', false); 
    
    threat_counters = horzcat(industry_threat_counts{:});
    
    if strcmp(threat_type, 'domestic')
        counters_check = build_domestic_counter(iucn_data_object, status_inds_to_use);
    elseif strcmp(threat_type, 'global')
        counters_check = build_global_counter(iucn_data_object, status_inds_to_use);
    end
    discrim = abs(sum(threat_counters - counters_check));
    if any(discrim > 1e-10)
        error(['mismatch in ', threat_type, ' satellite and records'])
    else 
       display(['all ', threat_type, 'records checked and correct']) 
    end
end


function [current_satellite, satellite_params] = sort_satellite_by_country(current_satellite, iucn_data_object, satellite_params) 
    
    industry_threat_counts = cellfun(@(x) sum(x, 2), current_satellite, 'UniformOutput', false); 

        if satellite_params.use_endemics == true
            sector_vec = build_sector_vec(satellite_params.sector_lengths);
            threat_matrix = horzcat(industry_threat_counts{:});
            threat_matrix = threat_matrix > 0;
            discrim = sum(threat_matrix, 2) == 1; 
            for country_ind = 1:numel(sector_vec)
                current_satellite{country_ind} = current_satellite{country_ind}.*repmat(discrim, [1, satellite_params.sector_lengths(country_ind)]);
            end
            industry_threat_counts = cellfun(@(x) sum(x, 2), current_satellite, 'UniformOutput', false); 
        end
     
        if strcmp(satellite_params.country_sort_type, 'sector_length') % OK %
            metric_to_sort = satellite_params.sector_lengths;
        elseif strcmp(satellite_params.country_sort_type, 'species_count') % OK %
            metric_to_sort = cellfun(@(x) sum( (x > 0), 1), industry_threat_counts);
        elseif strcmp(satellite_params.country_sort_type, 'threat_count') % OK %
            metric_to_sort = cellfun(@(x) sum(x(:)), industry_threat_counts);
        elseif strcmp(satellite_params.country_sort_type, 'threat_max')
            metric_to_sort = cellfun(@(x) max(x(:)), industry_threat_counts);
        elseif strcmp(satellite_params.country_sort_type, 'threat_mean')
            metric_to_sort = cellfun(@(x) mean(x(x>0)), industry_threat_counts);
        elseif strcmp(satellite_params.country_sort_type, 'GDP')
            metric_to_sort = cellfun(@(x) sum(x(:)), iucn_data_object.x); 
        end
        
        [~, sorted_country_inds] = sort(metric_to_sort);
 
        current_satellite = current_satellite(sorted_country_inds);
        satellite_params.sorted_country_names = iucn_data_object.iucn_country_names(sorted_country_inds);
        satellite_params.sector_lengths = sector_lengths(sorted_country_inds);
end

% (species, country, iucn_threat, iucn_status, industry_sector)
% (species, industry_sector, iucn_threat, iucn_status)

function [domestic_satellite, satellite_params] = build_domestic_satellite(iucn_data_object, satellite_params, species_characteristics)
    
    disp('building domestic satellite...')

    if (satellite_params.write_satellite_to_disk == true) && (~exist(satellite_params.satellite_filepath, 'dir'))
        mkdir([satellite_params.satellite_filepath]); 
    end

   domestic_satellite = aggregate_iucn_tensor('domestic', satellite_params.domestic_threats_to_aggregate, iucn_data_object, satellite_params, species_characteristics);
     
   [domestic_satellite, satellite_params] = sort_satellite(domestic_satellite, iucn_data_object, satellite_params);
        
   disp(['domestic satellite built at ', num2str(toc)])
          
end


function global_satellite = build_global_satellite(iucn_data_object, satellite_params, species_characteristics)
    
    disp('building global satellite...')
    global_threats = aggregate_iucn_tensor('global', satellite_params.global_threats_to_aggregate, iucn_data_object, satellite_params, species_characteristics);
    global_threats = sum([global_threats{:}], 2);
    species_num = satellite_params.classification_num;
    ghg_block = cellfun(@(x) (repmat(x', [species_num, 1])/iucn_data_object.global_ghg), iucn_data_object.ghg, 'UniformOutput', false);
    
    global_satellite = cellfun(@(x) x .* repmat(global_threats, [1, size(x, 2)]), ghg_block,'UniformOutput', false);
    
    [global_satellite, ~] = sort_satellite(global_satellite, iucn_data_object, satellite_params);
    
    disp(['global satellite built at ', num2str(toc)])
    
end

function current_satellite = aggregate_iucn_tensor(satellite_type, threats_to_aggregate, iucn_data_object, satellite_params, species_characteristics)
    
    if (satellite_params.return_satellite == true)
        
        if strcmp(satellite_type, 'domestic')
            current_satellite = arrayfun(@(x) zeros(satellite_params.classification_num, x), satellite_params.sector_lengths, 'UniformOutput', false);  
        else current_satellite = arrayfun(@(x) zeros(satellite_params.classification_num, 1), satellite_params.sector_lengths, 'UniformOutput', false); 
        end
        
    else
        disp(['writing ', satellite_type, ' satellite files to ', satellite_params.satellite_filepath])
        current_satellite = cell(1);
    end
    
     if strcmp(satellite_params.tensor_scale, 'country')
        country_indexes_to_use = 1:iucn_data_object.ncoun;
        
        for country_index = country_indexes_to_use
            
            load([satellite_params.pre_processed_tensor_file_prefix, iucn_data_object.iucn_country_code_names{country_index}, '.mat'])
       
            collapsed_tensor = collapse_species_tensor(satellite_params, species_characteristics, current_iucn_tensor, threats_to_aggregate);
        
            if strcmp(satellite_type, 'domestic')
                current_sector_length = satellite_params.sector_lengths(country_index);
            
                current_aggregated_country = zeros(satellite_params.classification_num, current_sector_length);
                current_aggregated_country(:, :) = collapsed_tensor(1, :, 1:current_sector_length);
                   
            else 
                current_aggregated_country = safe_collapse(collapsed_tensor, 1:size(collapsed_tensor, 2), 3); %collapse over industries
            end

            if (satellite_params.write_satellite_to_disk == true) 
               
                if strcmp(satellite_params.output_file_type, 'mat')
                    save([satellite_params.output_satellite_dir, satellite_params.system_type, '_', iucn_data_object.iucn_country_code_names{country_index}, '_',  'satellite.mat'], 'current_aggregated_country')
                elseif strcmp(satellite_params.output_file_type, 'csv')
                    csvwrite([satellite_params.output_satellite_dir, satellite_params.system_type, '_', iucn_data_object.iucn_country_code_names{country_index}, '_',  'satellite.csv'], current_aggregated_country)
                end
                
            end
        
            if satellite_params.return_satellite == true
                current_satellite{country_index}(:) = current_aggregated_country;
            end
            
            disp(['aggregated ' iucn_data_object.iucn_country_code_names{country_index} ' ' satellite_type ' data at ' num2str(toc)])
        
        end
        
     else   
        
        load([satellite_params.tensor_folder, 'iucn_tensor.mat'])   
        
        collapsed_tensor = collapse_species_tensor(satellite_params, species_characteristics, current_iucn_tensor, threats_to_aggregate);
        
        if strcmp(satellite_type, 'domestic')
            for country_ind = 1:iucn_data_object.ncoun
                current_satellite{country_ind}(:, :) = collapsed_tensor(country_ind, :, 1:satellite_params.sector_lengths(country_ind));
            end
        else 
            current_satellite = safe_collapse(current_iucn_tensor, 1:size(collapsed_tensor, 2), 3);      
        end
        
    end
        
end

function collapsed_tensor = collapse_species_tensor(satellite_params, species_characteristics, current_iucn_tensor, threats_to_aggregate)
        
        collapsed_tensor = safe_collapse(current_iucn_tensor(:, :, :, :, threats_to_aggregate), threats_to_aggregate, 5); % collapse over current iucn status levels   
        collapsed_tensor = safe_collapse(collapsed_tensor(:, :, :, satellite_params.status_inds_to_use), satellite_params.status_inds_to_use, 4); %collapse over current threats
        
        if numel(species_characteristics.sorted_q_row) > 1
            collapsed_tensor = collapsed_tensor(:, species_characteristics.sorted_q_row, :); % sort and select out species set
        else 
            collapsed_tensor_tmp = sptensor([size(collapsed_tensor, 1) 1 size(collapsed_tensor, 3)]);
            collapsed_tensor_tmp(:, 1, :) = collapsed_tensor(:, species_characteristics.sorted_q_row, :);
            collapsed_tensor = collapsed_tensor_tmp;
        end
    
        if satellite_params.collapse_through_species_sort_type == true
            collapsed_tensor = collapse_species(collapsed_tensor, satellite_params);
        end
      
end

function sector_list = build_sector_list(sector_lengths)
    sector_list = cell(size(sector_lengths));
    for country_ind = 1:length(sector_lengths);
        sector_list{country_ind} = country_ind*ones(sector_lengths(country_ind), 1);
    end
    sector_list = vertcat(sector_list{:});
end
% 
% function [sector_vec] = build_sector_vec(sector_lengths)
%     current_ncoun = length(sector_lengths);
%     sector_vec = zeros(1, current_ncoun);
%     for country_ind = 1:current_ncoun
%         sector_vec(country_ind) = sum(sector_lengths(1:(country_ind)));
%     end
% end
% 
% function [country_labels, sector_ticks] = build_country_labels(countries_to_label, sector_lengths, iucn_country_names)
%     sector_vec = build_sector_vec(sector_lengths);
%     if strcmp(countries_to_label, 'all')
%         sectors_to_label = find(sector_lengths > 0);
%     else 
%         sectors_to_label = find(ismember(iucn_country_names, countries_to_label));
%     end   
%     country_labels = iucn_country_names(sectors_to_label);
%     sector_vec_to_use = [0, sector_vec];
%     sector_ticks = floor(0.5*(sector_vec_to_use(sectors_to_label) + sector_vec_to_use(sectors_to_label + 1)));
% end
% 
% function [species_labels_to_use, species_ticks_to_use] = build_species_labels(classes_to_use, species_to_label)
%     [species_labels, species_class_counter] = unique(classes_to_use, 'stable');
%     if strcmp(species_to_label, 'all')
%         species_to_label = species_labels;
%     end
%     [~, ~, species_inds_to_use] = intersect(species_to_label, species_labels);
%     species_class_counter = vertcat(species_class_counter, numel(classes_to_use));
%     species_ticks_to_use = floor(0.5*(species_class_counter(species_inds_to_use) + species_class_counter(species_inds_to_use + 1)));
%     species_labels_to_use = species_labels(species_inds_to_use);
% end
% 
% 
% function display_satellite(satellite_to_use, countries_to_label, sector_lengths, iucn_country_names, classes_to_use, species_to_label)
% 
%     [country_labels, sector_ticks] = build_country_labels(countries_to_label, sector_lengths, iucn_country_names);
%     [species_labels_to_use, species_ticks_to_use] = build_species_labels(classes_to_use, species_to_label);
%     
%     imagesc(log(satellite_to_use)) 
%     set(gca, 'ytick', species_ticks_to_use)
%     set(gca,'YtickLabel', species_labels_to_use)
%     set(gca, 'xtick', sector_ticks)
%     set(gca, 'XtickLabel', country_labels)
%     set(gca, 'XTickLabelRotation', 90)
%     
% end




function [species_characteristics] = build_species_characteristics(iucn_data_object, satellite_params)
    display(['building satellite species characteristics'])
    if strcmp(satellite_params.species_taxons_to_use, 'all')
        species_taxons_to_use = unique(iucn_data_object.iucn_taxons_list);
    else species_taxons_to_use = satellite_params.species_taxons_to_use;
    end
    
    species_characteristics = struct();
    
    [~, ~, species_category_indexes] = intersect(species_taxons_to_use, iucn_data_object.iucn_threat_taxons);
    
    iucn_sort_category = select_class_system(iucn_data_object, satellite_params.species_sort_type, species_category_indexes);
    
    if (~strcmp(satellite_params.species_sub_category_to_use, 'all'))
        iucn_sub_category = select_class_system(iucn_data_object, satellite_params.species_sub_category_type, species_category_indexes);
        species_sub_category_to_use = find(ismember(iucn_sub_category, satellite_params.species_sub_category_to_use)); 
        [species_characteristics.sorted_species_classification, sorted_species_inds] = sort(iucn_sort_category(species_sub_category_to_use));       %sort remaining species by class type
        sorted_species_inds = species_sub_category_to_use(sorted_species_inds);
    else
        [species_characteristics.sorted_species_classification, sorted_species_inds] = sort(iucn_sort_category);       %sort remaining species by class type
    end
    
    sorted_species_indexes = species_category_indexes(sorted_species_inds);
    species_characteristics.sorted_species_taxons = iucn_data_object.iucn_threat_taxons(sorted_species_indexes);  
    species_characteristics.sorted_species_names = iucn_data_object.iucn_species_name(sorted_species_indexes); 
    
    [unique_species, q_row_indexes_to_use] = unique(iucn_data_object.iucn_taxons_list, 'stable');
    
    [~, ~, tmp_indexes] = intersect(species_characteristics.sorted_species_taxons, unique_species, 'stable');
    
    unique_q_row = iucn_data_object.q_row_list(q_row_indexes_to_use);
    species_characteristics.sorted_q_row = unique_q_row(tmp_indexes);
    
    species_characteristics.species_num = numel(species_characteristics.sorted_species_taxons);
        
end
  

function iucn_class = select_class_system(iucn_data_object, species_sort_type, inds_to_use)
    if strcmp(species_sort_type, 'species_kingdom') 
        iucn_class = iucn_data_object.iucn_species_kingdom(inds_to_use);
    elseif strcmp(species_sort_type, 'species_phylum') 
        iucn_class = iucn_data_object.iucn_species_phylum(inds_to_use);
    elseif strcmp(species_sort_type, 'species_class') 
        iucn_class = iucn_data_object.iucn_species_class(inds_to_use);
    elseif strcmp(species_sort_type, 'species_order') 
        iucn_class = iucn_data_object.iucn_species_order(inds_to_use);   
    elseif strcmp(species_sort_type, 'species_family') 
        iucn_class = iucn_data_object.iucn_species_family(inds_to_use);
    elseif strcmp(species_sort_type, 'species_genus') 
        iucn_class = iucn_data_object.iucn_species_genus(inds_to_use);
    elseif strcmp(species_sort_type, 'species_status') 
        iucn_class = iucn_data_object.iucn_threat_status(inds_to_use);
    end
end 


function [intCauseCounter] = build_global_counter(iucn_data_object, status_inds_to_use)
    elements_to_use = ~cellfun('isempty', iucn_data_object.iucn_vals);
    threats_to_use = ismember(iucn_data_object.threat_indexes_list, find(iucn_data_object.allcountriesflag)); %exclude global threats
    iucn_status_inds_to_use = ismember(iucn_data_object.iucn_status_indexes_list, status_inds_to_use);
    rows_to_use = find(elements_to_use & threats_to_use & iucn_status_inds_to_use)';
    intCauseCounter = build_counter(rows_to_use, species_characteristics.species_num, iucn_data_object.ncoun, iucn_data_object.q_row_list, iucn_data_object.country_indexes_list);
    
end

function [domCauseCounter] = build_domestic_counter(iucn_data_object, status_inds_to_use)
    elements_to_use = ~cellfun('isempty', iucn_data_object.iucn_vals);
    threats_to_use = ismember(iucn_data_object.threat_indexes_list, find(~iucn_data_object.allcountriesflag)); %exclude global threats
    iucn_status_inds_to_use = ismember(iucn_data_object.iucn_status_indexes_list, status_inds_to_use);
    rows_to_use = find(elements_to_use & threats_to_use & iucn_status_inds_to_use)';
    domCauseCounter = build_counter(rows_to_use, species_characteristics.species_num, iucn_data_object.ncoun, iucn_data_object.q_row_list, iucn_data_object.country_indexes_list);

end


function counter = build_counter(rows_to_use, nspec, ncoun, q_row_list, country_indexes_list)
    counter = zeros(nspec, ncoun);

    for row_ind = rows_to_use
        sp = q_row_list(row_ind);
    	cn = country_indexes_list(row_ind);
    	counter(sp,cn) = counter(sp,cn) + 1;
    end
end





