function [satellite_object] = build_iucn_satellite(satellite_inputs, satellite_params, system_type)
    
    display(['building satellite to ' system_type ' specification'])
       
    if (satellite_params.build_direct_satellite == true)    
        [satellite_object.direct_satellite, satellite_object.satellite_characteristics] = build_direct_satellite(satellite_params, satellite_inputs, system_type);
    end
    
    if (satellite_params.build_greenhouse_satellite == true)
        satellite_object.greenhouse_satellite = build_greenhouse_satellite(satellite_params, satellite_inputs, system_type);  
    end
    
    if ~strcmp(satellite_params.satellite_collapse_type, 'none')
        satellite_object.direct_satellite = run_collapse_satellite_routines('direct', satellite_params, satellite_inputs.satellite_characteristics, satellite_inputs.species_characteristics, system_type);
        satellite_object.greenhouse_satellite = run_collapse_satellite_routines('greenhouse', satellite_params, satellite_object.satellite_characteristics, satellite_inputs.species_characteristics, system_type);
   
    end
    
    if satellite_params.return_satellite_array == true
        satellite_object.direct_satellite = [satellite_object.direct_satellite{:}];
        satellite_object.greenhouse_satellite = [satellite_object.greenhouse_satellite{:}];
    end
    
    % display_satellite(satellite_object, satellite_object.satellite_characteristics, satellite_object.species_characteristics, ...
%                   satellite_object.satellite_characteristics.sector_lengths, satellite_object.satellite_characteristics.sorted_country_names)
    
    display(['satellite build finished at ', num2str(toc)])
    
end

function collapsed_satellite = run_collapse_satellite_routines(satellite_params, satellite_inputs, system_type, satellite_type)
    
    if strcmp(satellite_params.satellite_collapse_type, 'eora25')
        satellite_concordance = csvread([satellite_params.satellite_collapse_concordance_filename]);
        satellite_concordance = satellite_concordance';
        sector_classification_num = repmat(25, [numel(satellite_inputs.country_characteristics.country_codes), 1]);
    end
    
    if satellite_params.return_satellite == true
        collapsed_satellite = arrayfun(@(x) zeros(satellite_inputs.satellite_characteristics.classification_num, x), sector_classification_num, 'UniformOutput', false);  
    else collapsed_satellite = struct();
    end
    
    disp(['collapsing to ', satellite_params.satellite_collapse_type, '...'])
    
    for country_index = 1:numel(satellite_inputs.country_characteristics.country_codes)
        
        if strcmp(satellite_params.output_file_type, 'mat')
        	current_country = load([satellite_params.satellite_filepath, system_type, '_', satellite_inputs.country_characteristics.country_codes{country_index}, ...
                                    satellite_type '_satellite.mat']);
            fname = fieldnames(current_country);
            current_country = current_country.(fname{1});
            
        elseif strcmp(satellite_params.output_file_type, 'csv') 
            csvread([satellite_params.satellite_filepath, system_type, '_', satellite_inputs.country_characteristics.country_codes{country_index}, ...
                     satellite_type '_satellite.csv'])
        end
               
        if satellite_params.return_satellite == true
        	collapsed_satellite{country_index} = current_country*satellite_concordance;
        end
        
        disp([satellite_inputs.country_characteristics.country_codes{country_index}, ' collapsed at ' num2str(toc)])
        
    end
    
end


function collapsed_tensor = safe_collapse(current_tensor, collapse_inds, dim_to_collapse)

    if numel(collapse_inds) == 1
        collapsed_tensor = current_tensor;
    else collapsed_tensor = collapse(current_tensor, dim_to_collapse);
    end
    
end


function collapsed_species_tensor = collapse_species(iucn_tensor, satellite_params, satellite_characteristics)
    
    sz = size(iucn_tensor);
    unique_species_inds = satellite_params.unique_species_inds;
    collapsed_species_tensor = sptensor([sz(1) satellite_characteristics.classification_num sz(3:end)]);
    
    for collapse_ind = 1:satellite_characteristics.classification_num
        if collapse_ind < satellite_characteristics.classification_num
            species_set_to_use = unique_species_inds(collapse_ind):unique_species_inds(collapse_ind + 1);
            current_collapsed_tensor = safe_collapse(iucn_tensor(:, species_set_to_use, :), species_set_to_use, 2);
        else 
            species_set_to_use = unique_species_inds(collapse_ind):size(iucn_tensor,1);
            current_collapsed_tensor = safe_collapse(iucn_tensor(:, species_set_to_use, :), species_set_to_use, 2);
        end
        collapsed_species_tensor(:, collapse_ind, :) = current_collapsed_tensor;
    end  

end


function [current_satellite, satellite_characteristics] = sort_satellite(current_satellite, satellite_params, satellite_inputs)
                                                                                                                                                
    if strcmp(satellite_params.country_sort_type, 'none') 
        satellite_characteristics.sorted_country_names = satellite_inputs.country_characteristics.country_names;
    else
        if strcmp(satellite_params.country_sort_type, 'eora')
            disp('sorting to eora specification');
            [current_satellite, satellite_characteristics] = sort_satellite_to_eora(current_satellite, satellite_inputs.satellite_characteristics, satellite_inputs.country_characteristics.country_names, ...
                                                                                    satellite_inputs.satellite_characteristics.classification_num);
        else
            [current_satellite, satellite_characteristics] = sort_satellite_by_country(current_satellite, satellite_params, satellite_inputs.satellite_characteristics);
        end
    end
    
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

function [expanded_satellite, satellite_characteristics] = sort_satellite_to_eora(current_satellite, satellite_characteristics, iucn_country_names, classification_num) 
    
    eora_mapping_vector = satellite_characteristics.mapped_country_characteristics.mapping_vector;
    
    dummy_inds = find(eora_mapping_vector == 0); 
    dummy_country_names = satellite_characteristics.mapped_country_characteristics.country_names(dummy_inds);    
    expanded_satellite = cell(1, length(eora_mapping_vector));
    expanded_satellite((eora_mapping_vector > 0)) = current_satellite(eora_mapping_vector(eora_mapping_vector > 0));
    sector_lengths = cellfun(@(x) size(x, 2), expanded_satellite);
    sector_matches = (sector_lengths' == satellite_characteristics.mapped_country_characteristics.industry_specification(:, 1));
        
    dummy_satellite = build_dummy_satellite(eora_mapping_vector, sector_matches, classification_num, satellite_characteristics.mapped_country_characteristics.industry_specification);
    
    expanded_satellite = cellfun(@(x, y) [x y], expanded_satellite, dummy_satellite, 'UniformOutput', false);
    
    satellite_characteristics.sorted_country_names(dummy_inds) = dummy_country_names;
    satellite_characteristics.sorted_country_names(eora_mapping_vector > 0) = iucn_country_names(eora_mapping_vector(eora_mapping_vector > 0));

    satellite_characteristics.sector_lengths = cellfun(@(x) size(x, 2), expanded_satellite);
    
end
    



function [current_satellite, satellite_characteristics] = sort_satellite_by_country(current_satellite, satellite_params, satellite_characteristics) 
    
    industry_threat_counts = cellfun(@(x) sum(x, 2), current_satellite, 'UniformOutput', false); 

        if satellite_params.use_endemics == true
            sector_vec = build_sector_vec(satellite_characteristics.sector_lengths);
            threat_matrix = horzcat(industry_threat_counts{:});
            threat_matrix = threat_matrix > 0;
            discrim = sum(threat_matrix, 2) == 1; 
            for country_ind = 1:numel(sector_vec)
                current_satellite{country_ind} = current_satellite{country_ind}.*repmat(discrim, [1, satellite_characteristics.sector_lengths(country_ind)]);
            end
            industry_threat_counts = cellfun(@(x) sum(x, 2), current_satellite, 'UniformOutput', false); 
        end
     
        if strcmp(satellite_params.country_sort_type, 'sector_length') % OK %
            metric_to_sort = satellite_characteristics.sector_lengths;
        elseif strcmp(satellite_params.country_sort_type, 'species_count') % OK %
            metric_to_sort = cellfun(@(x) sum( (x > 0), 1), industry_threat_counts);
        elseif strcmp(satellite_params.country_sort_type, 'threat_count') % OK %
            metric_to_sort = cellfun(@(x) sum(x(:)), industry_threat_counts);
        elseif strcmp(satellite_params.country_sort_type, 'threat_max')
            metric_to_sort = cellfun(@(x) max(x(:)), industry_threat_counts);
        elseif strcmp(satellite_params.country_sort_type, 'threat_mean')
            metric_to_sort = cellfun(@(x) mean(x(x>0)), industry_threat_counts);
        elseif strcmp(satellite_params.country_sort_type, 'GDP')
           % metric_to_sort = cellfun(@(x) sum(x(:)), industry_characteristics.x); 
        end
        
        [~, sorted_country_inds] = sort(metric_to_sort);
 
        current_satellite = current_satellite(sorted_country_inds);
        satellite_characteristics.sorted_country_names = satellite_characteristics.country_names(sorted_country_inds);
        satellite_characteristics.sector_lengths = sector_lengths(sorted_country_inds);
end

% (species, country, iucn_threat, iucn_status, industry_sector)
% (species, industry_sector, iucn_threat, iucn_status)

function [direct_satellite, satellite_characteristics] = build_direct_satellite(satellite_params, satellite_inputs, system_type)
    
    disp('building direct satellite...')

    if (satellite_params.write_satellite_to_disk == true) && (~exist(satellite_params.satellite_filepath, 'dir'))
        mkdir([satellite_params.satellite_filepath]); 
    end

   direct_satellite = aggregate_iucn_tensor('direct', satellite_inputs.satellite_characteristics.direct_threats, satellite_params, satellite_inputs, system_type);
     
   [direct_satellite, satellite_characteristics] = sort_satellite(direct_satellite, satellite_params, satellite_inputs);
        
   disp(['direct satellite built at ', num2str(toc)])
          
end


function greenhouse_satellite = build_greenhouse_satellite(satellite_params, satellite_inputs, system_type)
    
    disp('building global satellite...')
    greenhouse_threats = aggregate_iucn_tensor('greenhouse', satellite_inputs.satellite_characteristics.greenhouse_threats, satellite_params, satellite_inputs, system_type);
    greenhouse_threats = sum([greenhouse_threats{:}], 2);
    species_num = satellite_inputs.satellite_characteristics.classification_num;
    ghg_block = cellfun(@(x) (repmat(x', [species_num, 1])/satellite_inputs.country_characteristics.global_ghg), satellite_inputs.country_characteristics.ghg, 'UniformOutput', false);
    
    greenhouse_satellite = cellfun(@(x) x .* repmat(greenhouse_threats, [1, size(x, 2)]), ghg_block,'UniformOutput', false);
    
    [greenhouse_satellite, ~] = sort_satellite(greenhouse_satellite, satellite_params, satellite_inputs);
    
    disp(['global satellite built at ', num2str(toc)])
    
end

function current_satellite = aggregate_iucn_tensor(satellite_type, threats_to_aggregate, satellite_params, satellite_inputs, system_type)
    
    if (satellite_params.return_satellite == true)
        
        if strcmp(satellite_type, 'direct')
            current_satellite = arrayfun(@(x) zeros(satellite_inputs.satellite_characteristics.classification_num, x), satellite_inputs.satellite_characteristics.sector_lengths, 'UniformOutput', false);  
        else
            current_satellite = arrayfun(@(x) zeros(satellite_inputs.satellite_characteristics.classification_num, 1), satellite_inputs.satellite_characteristics.sector_lengths, 'UniformOutput', false);
        end
        
    else
        disp(['writing ', satellite_type, ' satellite files to ', satellite_params.satellite_filepath])
        current_satellite = cell(1);
    end
    
     if strcmp(satellite_params.tensor_scale, 'country')
        
        country_indexes_to_use = 1:numel(satellite_inputs.country_characteristics.country_codes);
        
        for country_index = country_indexes_to_use
            
            load([satellite_params.iucn_tensor_file_prefix, satellite_inputs.country_characteristics.country_codes{country_index}, '.mat'])
            
            collapsed_tensor = collapse_species_tensor(satellite_params, satellite_inputs.satellite_characteristics, satellite_inputs.species_characteristics, current_iucn_tensor, threats_to_aggregate);
        
            if strcmp(satellite_type, 'direct')
                current_sector_length = satellite_inputs.satellite_characteristics.sector_lengths(country_index);
            
                current_aggregated_country = zeros(satellite_inputs.satellite_characteristics.classification_num, current_sector_length);
                current_aggregated_country(:, :) = collapsed_tensor(1, :, 1:current_sector_length);
                   
            else 
                current_aggregated_country = safe_collapse(collapsed_tensor, 1:size(collapsed_tensor, 2), 3); %collapse over industries
            end

            if (satellite_params.write_satellite_to_disk == true) 
               
                if strcmp(satellite_params.output_file_type, 'mat')
                    save([satellite_params.output_satellite_filepath, system_type, '_', satellite_inputs.country_characteristics.country_codes{country_index}, '_', satellite_type, '_satellite.mat'], 'current_aggregated_country')
                elseif strcmp(satellite_params.output_file_type, 'csv')
                    csvwrite([satellite_params.output_satellite_filepath, system_type, '_', satellite_inputs.country_characteristics.country_codes{country_index}, '_',  satellite_type, '_satellite.csv'], 'current_aggregated_country')
                end
                
            end
        
            if satellite_params.return_satellite == true
                current_satellite{country_index}(:) = current_aggregated_country;
            end
            
            disp([satellite_inputs.country_characteristics.country_codes{country_index} ' ' satellite_type ' threats satellite built at ' num2str(toc)])
        
        end
        
     else   
        
        load([satellite_params.tensor_folder, 'iucn_tensor.mat'])   
        
        collapsed_tensor = collapse_species_tensor(satellite_params, satellite_inputs.satellite_characteristics, satellite_inputs.species_characteristics, current_iucn_tensor, threats_to_aggregate);
        
        if strcmp(satellite_type, 'direct')
            for country_ind = 1:numel(satellite_inputs.country_characteristics.country_codes);
                current_satellite{country_ind}(:, :) = collapsed_tensor(country_ind, :, 1:satellite_inputs.satellite_characteristics.sector_lengths(country_ind));
            end
        else 
            current_satellite = safe_collapse(current_iucn_tensor, 1:size(collapsed_tensor, 2), 3);      
        end
        
    end
        
end

function collapsed_tensor = collapse_species_tensor(satellite_params, satellite_characteristics, species_characteristics, current_iucn_tensor, threats_to_aggregate)
        
        collapsed_tensor = safe_collapse(current_iucn_tensor(:, :, :, :, threats_to_aggregate), threats_to_aggregate, 5); % collapse over current iucn status levels   
        collapsed_tensor = safe_collapse(collapsed_tensor(:, :, :, satellite_characteristics.status_inds_to_use), satellite_characteristics.status_inds_to_use, 4); %collapse over current threats
        
        if numel(species_characteristics.sorted_q_row) > 1
            collapsed_tensor = collapsed_tensor(:, species_characteristics.sorted_q_row, :); % sort and select out species set
        else 
            collapsed_tensor_tmp = sptensor([size(collapsed_tensor, 1) 1 size(collapsed_tensor, 3)]);
            collapsed_tensor_tmp(:, 1, :) = collapsed_tensor(:, species_characteristics.sorted_q_row, :);
            collapsed_tensor = collapsed_tensor_tmp;
        end
    
        if satellite_params.collapse_through_species_sort_type == true
            collapsed_tensor = collapse_species(collapsed_tensor, satellite_characteristics, satellite_params);
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



%%%%%%%    SATELLITE CHECKING ROUTINES 

% function check_satellite(split_satellite, processed_iucn_data, status_inds_to_use, threat_type)
%     industry_threat_counts = cellfun(@(x) sum(x, 2), split_satellite, 'UniformOutput', false); 
%     
%     threat_counters = horzcat(industry_threat_counts{:});
%     
%     if strcmp(threat_type, 'direct')
%         counters_check = build_direct_counter(processed_iucn_data, status_inds_to_use);
%     elseif strcmp(threat_type, 'greenhouse')
%         counters_check = build_greenhouse_counter(processed_iucn_data, status_inds_to_use);
%     end
%     discrim = abs(sum(threat_counters - counters_check));
%     if any(discrim > 1e-10)
%         error(['mismatch in ', threat_type, ' satellite and records'])
%     else 
%        display(['all ', threat_type, 'records checked and correct']) 
%     end
% end
% function [intCauseCounter] = build_greenhouse_counter(processed_iucn_data, status_inds_to_use)
%     elements_to_use = ~cellfun('isempty', processed_iucn_data.iucn_vals);
%     threats_to_use = ismember(processed_iucn_data.threat_indexes_list, find(processed_iucn_data.allcountriesflag)); %exclude global threats
%     iucn_status_inds_to_use = ismember(processed_iucn_data.iucn_status_indexes_list, status_inds_to_use);
%     rows_to_use = find(elements_to_use & threats_to_use & iucn_status_inds_to_use)';
%     intCauseCounter = build_counter(rows_to_use, species_characteristics.species_num, processed_iucn_data.ncoun, processed_iucn_data.q_row_list, processed_iucn_data.country_indexes_list);
%     
% end
% 
% function [domCauseCounter] = build_direct_counter(processed_iucn_data, status_inds_to_use)
%     elements_to_use = ~cellfun('isempty', processed_iucn_data.iucn_vals);
%     threats_to_use = ismember(processed_iucn_data.threat_indexes_list, find(~processed_iucn_data.allcountriesflag)); %exclude global threats
%     iucn_status_inds_to_use = ismember(processed_iucn_data.iucn_status_indexes_list, status_inds_to_use);
%     rows_to_use = find(elements_to_use & threats_to_use & iucn_status_inds_to_use)';
%     domCauseCounter = build_counter(rows_to_use, species_characteristics.species_num, processed_iucn_data.ncoun, processed_iucn_data.q_row_list, processed_iucn_data.country_indexes_list);
% 
% end

% 
% function counter = build_counter(rows_to_use, nspec, ncoun, q_row_list, country_indexes_list)
%     counter = zeros(nspec, ncoun);
% 
%     for row_ind = rows_to_use
%         sp = q_row_list(row_ind);
%     	cn = country_indexes_list(row_ind);
%     	counter(sp,cn) = counter(sp,cn) + 1;
%     end
% end





