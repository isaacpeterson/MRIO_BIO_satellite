function satellite_inputs = build_processed_iucn_data(industry_inputs, iucn_data_params, global_params, satellite_params)
    
    tic
    
    if (iucn_data_params.build_processed_iucn_data == false) && exist([iucn_data_params.processed_iucn_data_filename_prefix, global_params.system_type, '.mat'], 'file')
        disp(['loading processed iucn data from ', [iucn_data_params.processed_iucn_data_filename_prefix, global_params.system_type, '.mat']])
        processed_data = load([iucn_data_params.processed_iucn_data_filename_prefix, global_params.system_type, '.mat']);
        obj_names = fieldnames(processed_data);
        processed_data = processed_data.(obj_names{1});
    
    else 
    
        processed_data = struct();

        processed_data.iucn_data = build_iucn_data(iucn_data_params, industry_inputs.un_to_iucn_codes, industry_inputs.industry_characteristics, global_params);
    
        processed_data.industry_data = build_industry_data(processed_data.iucn_data, ...
                                                           iucn_data_params, ...
                                                           industry_inputs.un_to_iucn_codes, ...
                                                           industry_inputs.industry_characteristics, ...
                                                           global_params.system_type);
    
        if ~exist(global_params.processed_data_filepath, 'dir')
            mkdir(global_params.processed_data_filepath); 
        end
    
        if (iucn_data_params.save_processed_data == true)
            disp(['saving processed iucn data to ', global_params.processed_data_filepath])
            save([iucn_data_params.processed_iucn_data_filename_prefix, global_params.system_type, '.mat'], 'processed_data', '-v7.3');
        end
        
    end
   
    
    if (iucn_data_params.overwrite_tensors == true)
        disp(['processing and saving ', iucn_data_params.tensor_scale, ' level tensors...'])
        build_iucn_tensors(processed_data, iucn_data_params, global_params.system_type, industry_inputs.industry_characteristics, industry_inputs.un_to_iucn_codes);
    end
    
    disp('building satellite inputs ...')
    satellite_inputs = build_satellite_inputs(processed_data, satellite_params, industry_inputs, global_params.system_type);
    
end

function satellite_inputs = build_satellite_inputs(processed_iucn_data, satellite_params, industry_inputs, system_type)
    
    satellite_inputs = struct();
    
    satellite_inputs.species_characteristics = build_species_characteristics(processed_iucn_data, satellite_params);
    
    satellite_inputs.satellite_characteristics = build_satellite_parameters(satellite_params, processed_iucn_data, satellite_inputs.species_characteristics, ...
                                                                                industry_inputs.industry_characteristics, industry_inputs.un_to_iucn_codes);
    
    satellite_inputs.iucn_industry_characteristics = buld_iucn_industry_characteristics(system_type, industry_inputs.industry_characteristics, ...
                                                                                        processed_iucn_data.iucn_data.country_characteristics.country_codes, ...
                                                                                        industry_inputs.un_to_iucn_codes.iucn_country_codes, satellite_params.include_ghg); 
                                                                                    
    satellite_inputs.iucn_data_country_characteristics = processed_iucn_data.iucn_data.country_characteristics;
    
end 


function iucn_industry_characteristics = buld_iucn_industry_characteristics(system_type, industry_characteristics, country_codes, iucn_codes, include_ghg)

    iucn_industry_characteristics = struct();
    
    if strcmp(system_type, 'eora') 
        iucn_industry_characteristics.x_names = reorder_to_iucn(industry_characteristics.x_un_names, country_codes, iucn_codes);
    elseif strcmp(system_type, 'hscpc')   
        display(['fix this routine'])
    end
    
    if (include_ghg == true)
            
        if strcmp(system_type, 'eora') 
    
                iucn_industry_characteristics.ghg = reorder_to_iucn(industry_characteristics.ghg, country_codes, iucn_codes);
                iucn_industry_characteristics.global_ghg = sum(cell2mat(iucn_industry_characteristics.ghg)); 
                
        else 
            display(['HSCPC ghg data unavailable'])
        end
            
    end
    
end

function satellite_characteristics = build_satellite_parameters(satellite_params, processed_iucn_data, species_characteristics, industry_characteristics, un_to_iucn_codes)
    
    satellite_characteristics = struct();
    
    if strcmp(satellite_params.direct_threats, 'all')
        satellite_characteristics.direct_threats = find(~processed_iucn_data.iucn_data.threat_characteristics.greenhouse_threats)';
    else satellite_characteristics.direct_threats = satellite_params.direct_threats;
    end
    
    if strcmp(satellite_params.greenhouse_threats, 'all')
        satellite_characteristics.greenhouse_threats = find(processed_iucn_data.iucn_data.threat_characteristics.greenhouse_threats)';
    else satellite_characteristics.greenhouse_threats = satellite_params.greenhouse_threats;
    end   
    
    if  strcmp(satellite_params.status_levels_to_use, 'all')
       satellite_characteristics.status_inds_to_use = 1:length(processed_iucn_data.iucn_data.iucn_species_data.iucn_status_names);
    else
        [~, satellite_characteristics.status_inds_to_use] = intersect(processed_iucn_data.iucn_data.iucn_species_data.iucn_status_names, satellite_params.status_levels_to_use);
    end

    if satellite_params.collapse_through_species_sort_type == true
        [~, unique_species_inds] = unique(species_characteristics.sorted_species_classification, 'stable');
        satellite_characteristics.classification_num = numel(unique_species_inds);
    else satellite_characteristics.classification_num = numel(species_characteristics.sorted_q_row);
    end
    
    satellite_characteristics.sector_lengths = cellfun('length', processed_iucn_data.industry_data.x);
    
    if strcmp(satellite_params.country_sort_type, 'eora')
    	satellite_characteristics.mapped_country_characteristics = find_eora_industry_specification(industry_characteristics, un_to_iucn_codes, processed_iucn_data.iucn_data.country_characteristics.country_codes); 
    end  
    
end

function [satellite_country_characteristics] = find_eora_industry_specification(industry_characteristics, un_to_iucn_codes, iucn_country_code_names)
    
    satellite_country_characteristics = struct();
    
    [~, un_inds, iucn_inds] = intersect(un_to_iucn_codes.iucn_country_codes, iucn_country_code_names, 'stable');

    un_to_iucn_inds = zeros(length(un_to_iucn_codes.un_country_codes), 1);
    un_to_iucn_inds(un_inds) = iucn_inds;

    [~, ~, un_inds_to_use] = intersect(industry_characteristics.unique_country_codes, un_to_iucn_codes.un_country_codes, 'stable');
    
    satellite_country_characteristics.mapping_vector = un_to_iucn_inds(un_inds_to_use);

    satellite_country_characteristics.country_names = industry_characteristics.unique_countries;

    satellite_country_characteristics.industry_specification = zeros(numel(satellite_country_characteristics.mapping_vector), 2);
    
    for country_ind = 1:numel(satellite_country_characteristics.mapping_vector)
        current_country_data = strcmp(industry_characteristics.country_codes_list, industry_characteristics.unique_country_codes(country_ind));   
        industries_list = find(current_country_data & industry_characteristics.industry_data_list);
        commodities_list = find(current_country_data & industry_characteristics.commodity_data_list);
        satellite_country_characteristics.industry_specification(country_ind, :) = [numel(industries_list), numel(commodities_list)];
    end
    
  
end
    

function [species_characteristics] = build_species_characteristics(processed_iucn_data, satellite_params)
    
    display(['building satellite species characteristics'])
    
    if strcmp(satellite_params.species_taxons_to_use, 'all')
        species_taxons_to_use = unique(processed_iucn_data.iucn_data.iucn_species_data.iucn_taxons_list);
    else species_taxons_to_use = satellite_params.species_taxons_to_use;
    end
    
    [~, ~, species_category_indexes] = intersect(species_taxons_to_use, processed_iucn_data.iucn_data.iucn_threat_data.iucn_threat_taxons);
    
    iucn_sort_category = select_class_system(processed_iucn_data, satellite_params.species_sort_type, species_category_indexes);
    
    if (~strcmp(satellite_params.species_sub_category_to_use, 'all'))
        iucn_sub_category = select_class_system(processed_iucn_data, satellite_params.species_sub_category_type, species_category_indexes);
        species_sub_category_to_use = find(ismember(iucn_sub_category, satellite_params.species_sub_category_to_use)); 
        [sorted_species_classification, sorted_species_inds] = sort(iucn_sort_category(species_sub_category_to_use));       %sort remaining species by class type
        sorted_species_inds = species_sub_category_to_use(sorted_species_inds);
    else
        [sorted_species_classification, sorted_species_inds] = sort(iucn_sort_category);       %sort remaining species by class type
    end
    
    sorted_species_indexes = species_category_indexes(sorted_species_inds);
    
    species_characteristics = structfun(@(x) x(sorted_species_indexes), processed_iucn_data.iucn_data.iucn_threat_data, 'un', false);
     
    species_characteristics.sorted_species_classification = sorted_species_classification;
    
    [unique_species, q_row_indexes_to_use] = unique(processed_iucn_data.iucn_data.iucn_species_data.iucn_taxons_list, 'stable');
    
    [~, ~, tmp_indexes] = intersect(species_characteristics.iucn_threat_taxons, unique_species, 'stable');
    
    unique_q_row = processed_iucn_data.iucn_data.iucn_species_data.q_row_list(q_row_indexes_to_use);
    species_characteristics.sorted_q_row = unique_q_row(tmp_indexes);
    
    species_characteristics.species_num = numel(species_characteristics.iucn_threat_taxons);
        
end
  

function iucn_class = select_class_system(processed_iucn_data, species_sort_type, inds_to_use)
    if strcmp(species_sort_type, 'species_kingdom') 
        iucn_class = processed_iucn_data.iucn_data.iucn_threat_data.iucn_species_kingdom(inds_to_use);
    elseif strcmp(species_sort_type, 'species_phylum') 
        iucn_class = processed_iucn_data.iucn_data.iucn_threat_data.iucn_species_phylum(inds_to_use);
    elseif strcmp(species_sort_type, 'species_class') 
        iucn_class = processed_iucn_data.iucn_data.iucn_threat_data.iucn_species_class(inds_to_use);
    elseif strcmp(species_sort_type, 'species_order') 
        iucn_class = processed_iucn_data.iucn_data.iucn_threat_data.iucn_species_order(inds_to_use);   
    elseif strcmp(species_sort_type, 'species_family') 
        iucn_class = processed_iucn_data.iucn_data.iucn_threat_data.iucn_species_family(inds_to_use);
    elseif strcmp(species_sort_type, 'species_genus') 
        iucn_class = processed_iucn_data.iucn_data.iucn_threat_data.iucn_species_genus(inds_to_use);
    elseif strcmp(species_sort_type, 'species_status') 
        iucn_class = processed_iucn_data.iucn_data.iucn_threat_data.iucn_threat_status(inds_to_use);
    end
end 



function processed_industry_data = build_industry_data(processed_iucn_data, iucn_data_params, un_to_iucn_codes, industry_characteristics, system_type)
     
    if strcmp(iucn_data_params.iucn_data_type, 'new')
        iucn_codes = un_to_iucn_codes.iucn_country_codes;
    else    
        iucn_codes = un_to_iucn_codes.un_country_codes;
    end
    
    threat_concordance = build_threat_concordance(processed_iucn_data, iucn_data_params, un_to_iucn_codes, system_type);
    
    processed_industry_data.x = reorder_to_iucn(industry_characteristics.x_un, processed_iucn_data.country_characteristics.country_codes, iucn_codes);    
    
    if strcmp(system_type, 'eora')
       processed_industry_data.industry_indexes = build_eora_industry_indexes(processed_iucn_data.iucn_species_data.iucn_data_length, processed_iucn_data.iucn_species_data.threat_indexes_list, ...
                                                      processed_iucn_data.iucn_species_data.country_indexes_list, threat_concordance, processed_iucn_data.country_characteristics.ncoun);
    elseif strcmp(system_type, 'hscpc')
       processed_industry_data.industry_indexes = build_hscpc_industry_indexes(processed_iucn_data.iucn_species_data.iucn_data_length, processed_iucn_data.iucn_species_data.threat_indexes_list, ...
                                                       threat_concordance, processed_iucn_data.threat_characteristics.threat_num);
    end

    processed_industry_data.sector_lengths = cellfun('length', processed_industry_data.industry_indexes);
    processed_industry_data.scaled_industry_vals = scale_industry_vals(processed_iucn_data.iucn_species_data.iucn_data_length, ...
                                                                       processed_industry_data.x, ...
                                                                       processed_iucn_data.iucn_species_data.country_indexes_list, ...
                                                                       processed_industry_data.industry_indexes);  
                                                                     
end

function processed_iucn_data = build_iucn_data(iucn_data_params, un_to_iucn_codes, industry_characteristics, global_params)
    
    disp(['processing iucn data to ', global_params.system_type, ' specification...'])
        
    processed_iucn_data = struct();
        
   	if strcmp(iucn_data_params.iucn_data_type, 'new')  
        
        processed_iucn_data.iucn_threat_data = read_raw_iucn_threat_data(iucn_data_params.iucn_data_threats_filename);
    	processed_iucn_data.iucn_species_data = process_raw_iucn_species_data(iucn_data_params.iucn_data_species_filename, processed_iucn_data.iucn_threat_data, ...
                                                                              industry_characteristics.industry_codes_to_use, un_to_iucn_codes.iucn_country_codes);
    else
    	processed_iucn_data = initialise_object_from_old_data(iucn_data_params.old_iucn_data_filename);
    end

    disp(['processed raw iucn data at ', num2str(toc)])
    
    processed_iucn_data.country_characteristics = build_iucn_country_characteristics(iucn_data_params.read_iucn_countries_from_file, iucn_data_params.iucn_data_type, un_to_iucn_codes, ...
                                                                                     processed_iucn_data.iucn_species_data.iucn_country_codes_list, processed_iucn_data.iucn_species_data.country_names_list);
                                                                                 
    processed_iucn_data.threat_characteristics = build_threat_categories(iucn_data_params.old_threat_cause_class_filename, iucn_data_params.greenhouse_flag_filename, ...
                                                                iucn_data_params.read_threat_classification_from_file, processed_iucn_data.iucn_threat_data);
    
    processed_iucn_data = build_iucn_species_indexes(processed_iucn_data);
    
end

function threat_concordance = build_threat_concordance(processed_iucn_data, iucn_data_params, un_to_iucn_codes, system_type)
    
    disp('building threat concordances ...')
    
    if strcmp(system_type, 'hscpc')
        threat_concordance = build_hscpc_threat_concordance(iucn_data_params.hscpc_concordance_filename, processed_iucn_data.threat_characteristics, iucn_data_params.hscpc_sector_num);
                                                                                      
    elseif strcmp(system_type, 'eora') 
        threat_concordance = build_eora_threat_concordance(iucn_data_params.eora_concordance_file_prefix, processed_iucn_data.country_characteristics.country_codes, ...
                                                                                    un_to_iucn_codes, processed_iucn_data.country_characteristics.ncoun);
    end
end


function build_iucn_tensors(processed_data, iucn_data_params, system_type, industry_characteristics, un_to_iucn_codes)

    disp(['writing tensors to ' iucn_data_params.iucn_tensor_file_prefix])
    
    if ~exist(iucn_data_params.iucn_tensor_file_prefix, 'dir')
        mkdir(iucn_data_params.iucn_tensor_file_prefix); 
    end
    
    if strcmp(iucn_data_params.tensor_scale, 'country')
        
        for country_index = 1:processed_data.iucn_data.country_characteristics.ncoun
         
            rows_to_use = find(ismember(processed_data.iucn_data.iucn_species_data.country_indexes_list, country_index));
            
            current_iucn_tensor = build_current_tensor(processed_data, iucn_data_params.tensor_scale, rows_to_use);
            current_tensor_filename = [iucn_data_params.iucn_tensor_file_prefix, system_type, '_', processed_data.iucn_data.country_characteristics.country_codes{country_index}, '.mat'];
            save(current_tensor_filename, 'current_iucn_tensor')
            disp([current_tensor_filename, ' done at ' num2str(toc)])
            
        end   
        
    else
        
        current_iucn_tensor = build_current_tensor(processed_data, iucn_data_params.tensor_scale, 1:length(processed_data.country_indexes_list));
        disp(['global iucn tensor built at ' num2str(toc)])
        	
        current_tensor_filename = [iucn_data_params.iucn_tensor_file_prefix, 'global_iucn_tensor.mat'];
        save(current_tensor_filename, 'current_iucn_tensor')
              
    end
    
end

function current_iucn_tensor = build_current_tensor(processed_data, tensor_type, rows_to_use)

    current_tensor_block = zeros(sum(processed_data.industry_data.sector_lengths(rows_to_use)), 1);
    
    if strcmp(tensor_type, 'global')
        current_tensor_block(:, 1) = build_iucn_tensor_indexes(processed_data.country_indexes_list(rows_to_use), processed_data.industry_data.sector_lengths(rows_to_use));
    else current_tensor_block(:, 1) = ones(sum(processed_data.industry_data.sector_lengths(rows_to_use)), 1);
    end
        
    current_tensor_block(:, 2) = build_iucn_tensor_indexes(processed_data.iucn_data.iucn_species_data.q_row_list(rows_to_use), processed_data.industry_data.sector_lengths(rows_to_use));
    current_tensor_block(:, 3) = cat(2, processed_data.industry_data.industry_indexes{rows_to_use});
    current_tensor_block(:, 4) = build_iucn_tensor_indexes(processed_data.iucn_data.iucn_species_data.iucn_status_indexes_list(rows_to_use), processed_data.industry_data.sector_lengths(rows_to_use));
    current_tensor_block(:, 5) = build_iucn_tensor_indexes(processed_data.iucn_data.iucn_species_data.threat_indexes_list(rows_to_use), processed_data.industry_data.sector_lengths(rows_to_use));

    current_scaled_industry_vals = cat(1, processed_data.industry_data.scaled_industry_vals{rows_to_use});
    current_iucn_tensor = sptensor(double(current_tensor_block), double(current_scaled_industry_vals), double(max(current_tensor_block)));
   
end

function tensor_index_array = build_iucn_tensor_indexes(current_iucn_data, current_lengths)

    tensor_index_array = arrayfun(@(x, y) repmat(x, [y, 1]), current_iucn_data, current_lengths, 'UniformOutput', false);
    tensor_index_array = cell2mat(tensor_index_array); 
    
end


function iucn_country_characteristics = build_iucn_country_characteristics(read_iucn_countries_from_file, iucn_data_type, un_to_iucn_codes, country_codes_list, country_names_list)
    
    if read_iucn_countries_from_file
       
       if strcmp(iucn_data_type, 'old') 
            
            eora_codes = load_eora_country_codes(iucn_data_params.eora_countries_filename);
            country_names = eora_codes.eora_country_names;
            country_codes = eora_codes.eora_country_codes;
            empty_names = strcmp(country_codes, ' ');
           
       elseif strcmp(iucn_data_type, 'new')
           
            [country_codes, sort_inds] = sort(un_to_iucn_codes.iucn_country_codes);
            country_names = un_to_iucn_codes.country_names(sort_inds);
            empty_names = cellfun('isempty', country_codes);
       
       end
       
    else 
       
        [country_codes, unique_inds, ~] = unique(country_codes_list);
        country_names = country_names_list(unique_inds); 
        empty_names = cellfun('isempty', country_codes);

    end
    
    iucn_country_characteristics = struct();
    iucn_country_characteristics.country_codes = country_codes(~empty_names);
    iucn_country_characteristics.country_names = country_names(~empty_names);
    iucn_country_characteristics.ncoun = length(iucn_country_characteristics.country_codes);
    
end


function iucn_threat_characteristics = build_threat_categories(old_threat_cause_class_filename, greenhouse_flag_filename, read_threat_classification_from_file, iucn_threat_data)
    
    disp(['building threat classifications ...'])
    
    fid = fopen(old_threat_cause_class_filename);
        old_threat_cause_class_data = textscan(fid,'%s %s %f', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);

    iucn_threat_characteristics.old_threat_cause_class = strrep(old_threat_cause_class_data{1}, '''','');
    old_threat_cause_names = old_threat_cause_class_data{2};
    old_threat_cause_group = old_threat_cause_class_data{3};
        
    if read_threat_classification_from_file
        iucn_threat_characteristics.threat_cause_class = iucn_threat_characteristics.old_threat_cause_class;
        iucn_threat_characteristics.threat_cause_group = old_threat_cause_group;
   else
       iucn_threat_characteristics.threat_cause_class = unique(iucn_threat_data.iucn_threats);

       [~, name_inds, old_name_inds] = intersect(iucn_threat_characteristics.threat_cause_class, iucn_threat_characteristics.old_threat_cause_class, 'stable');
       threat_cause_names = cell(length(iucn_threat_characteristics.threat_cause_class), 1);
       threat_cause_names(name_inds) = old_threat_cause_names(old_name_inds);
       unlisted_set = setdiff(1:length(iucn_threat_characteristics.threat_cause_class), name_inds);
       threat_cause_names(unlisted_set) = iucn_threat_characteristics.threat_cause_class(unlisted_set);
       
       iucn_threat_characteristics.threat_cause_names = threat_cause_names;
       iucn_threat_characteristics.threat_cause_group = generate_threat_cause_group(iucn_threat_characteristics.threat_cause_class);
    end
    
    iucn_threat_characteristics.threat_num = length(iucn_threat_characteristics.threat_cause_class);
    
    greenhouse_flag = load(greenhouse_flag_filename);
    obj_name = fieldnames(greenhouse_flag);
    greenhouse_flag = greenhouse_flag.(obj_name{1});
    
    [~, old_threat_inds, new_threat_inds] = intersect(iucn_threat_characteristics.old_threat_cause_class, iucn_threat_characteristics.threat_cause_class);
    iucn_threat_characteristics.greenhouse_threats = zeros(size(iucn_threat_characteristics.threat_cause_class));
    iucn_threat_characteristics.greenhouse_threats(new_threat_inds) = greenhouse_flag(old_threat_inds);
    
end


function [eora_codes] = load_eora_country_codes(eora_countries_filename)
    
    [~,~, eora_country_data] = xlsread(eora_countries_filename); 
    eora_country_names = eora_country_data(3:end,3);
    eora_country_codes = eora_country_data(3:end,4);
    eora_country_codes = strrep(eora_country_codes, 'TZM', 'TZA');   %BIG NOTE - THIS WAS THE OTHER WAY AROUND IN THE PREVIOUS CODE - i.e. TANZANIA WAS EXCLUDED FROM THE ANALYSIS
    empty_names = strcmp(eora_country_codes, ' ');
    eora_country_codes = eora_country_codes(~empty_names);
    eora_country_names = eora_country_names(~empty_names);
    [eora_codes.eora_country_codes, sorted_inds] = sort(eora_country_codes);
    eora_codes.eora_country_names = eora_country_names(sorted_inds);
    
end   

    
function iucn_threat_data = read_raw_iucn_threat_data(new_iucn_data_threats_filename)

    fid = fopen(new_iucn_data_threats_filename);
        raw_iucn_threat_data = textscan(fid,'%f %f %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);

    empty_threats = cellfun('isempty', raw_iucn_threat_data{14}); 
    
    iucn_threat_data.iucn_threats = raw_iucn_threat_data{14}(~empty_threats);
    iucn_threat_data.iucn_threat_taxons = raw_iucn_threat_data{2}(~empty_threats);
    iucn_threat_data.iucn_species_kingdom = raw_iucn_threat_data{3}(~empty_threats);
    iucn_threat_data.iucn_species_phylum = raw_iucn_threat_data{4}(~empty_threats);
    iucn_threat_data.iucn_species_class = raw_iucn_threat_data{5}(~empty_threats);
    iucn_threat_data.iucn_species_order = raw_iucn_threat_data{6}(~empty_threats);
    iucn_threat_data.iucn_species_family = raw_iucn_threat_data{7}(~empty_threats);
    iucn_threat_data.iucn_species_genus = raw_iucn_threat_data{8}(~empty_threats);
    iucn_threat_data.iucn_species_name = raw_iucn_threat_data{9}(~empty_threats);
    iucn_threat_data.iucn_threat_status = raw_iucn_threat_data{10}(~empty_threats);
    
end
    
    
% process_raw_iucn_species_data(iucn_data_params.iucn_data_species_filename, processed_iucn_data.iucn_threat_data, ...
%                                                                               industry_characteristics.industry_codes_to_use, un_to_iucn_codes.iucn_country_codes);
                                                                          
function iucn_species_data = process_raw_iucn_species_data(new_iucn_data_species_filename, iucn_threat_data, industry_codes_to_use, iucn_country_codes)
    
    display('loading raw iucn data')
    fid = fopen(new_iucn_data_species_filename);
        data_block = textscan(fid,'%f %f %s %s %s %s %s %s %s %s %s %s %s %s %s', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);
    
    raw_species_data.taxons = data_block{2};
    raw_species_data.country_codes = data_block{12};
    raw_species_data.country_names = data_block{13};

    iucn_codes_to_use = unique(iucn_country_codes(industry_codes_to_use));  % remove codes from list that have no industry data
    iucn_codes_to_use = iucn_codes_to_use(~cellfun('isempty', iucn_codes_to_use));  %remove empty code names from list                                                                          
    iucn_species_list_inds_to_use = ismember(raw_species_data.country_codes, iucn_codes_to_use); %find all intersecting countries in iucn_list and UN list, i.e. remove iucn countries that there is no industrial data for.

    excluded_countries = unique(raw_species_data.country_names(~iucn_species_list_inds_to_use));
    
    raw_species_data = structfun(@(x) x(iucn_species_list_inds_to_use), raw_species_data, 'un', false);
    
    [~, sorted_iucn_code_inds] = sort(raw_species_data.country_codes);
    
    raw_species_data = structfun(@(x) x(sorted_iucn_code_inds), raw_species_data, 'un', false);
    
    row_length = length(raw_species_data.country_codes);

    display('processing raw iucn data')
    [expanded_iucn_threats, iucn_lengths] = build_iucn_threat_data(iucn_threat_data.iucn_threats, iucn_threat_data.iucn_threat_taxons, raw_species_data.taxons, row_length);    
    [expanded_iucn_status, ~] = build_iucn_threat_data(iucn_threat_data.iucn_threat_status, iucn_threat_data.iucn_threat_taxons, raw_species_data.taxons, row_length);
    expanded_iucn_status = strrep(expanded_iucn_status, '/', '_');
  
    iucn_species_data = struct();
    iucn_species_data.iucn_threats_list = expanded_iucn_threats;
    iucn_species_data.iucn_status_list = expanded_iucn_status;
    iucn_species_data.excluded_countries = excluded_countries;
    iucn_species_data.iucn_status_names = unique(iucn_species_data.iucn_status_list); % 11 extinction levels (VUlnerable,EXtinct, etc)    
    iucn_species_data.country_names_list = expand_iucn_species_data(raw_species_data.country_names, row_length, iucn_lengths);
    iucn_species_data.iucn_taxons_list = expand_iucn_species_data(raw_species_data.taxons, row_length, iucn_lengths);
    iucn_species_data.iucn_country_codes_list = expand_iucn_species_data(raw_species_data.country_codes, row_length, iucn_lengths);
    iucn_species_data.q_row_list = build_q_row(iucn_species_data.iucn_taxons_list);
    iucn_species_data.iucn_data_length = length(iucn_species_data.q_row_list);
  
end


function processed_iucn_data = initialise_object_from_old_data(old_iucn_data_filename)
    fid = fopen(old_iucn_data_filename);
        iucn_data = textscan(fid,'%f %s %s %s %s %s %s', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);

    processed_iucn_data.q_row_list = iucn_data{1};
    processed_iucn_data.country_names_list = iucn_data{2};
    processed_iucn_data.iucn_status_list = iucn_data{3}; % 11 threat levels (VUlnerable,EXtinct, etc)
    processed_iucn_data.iucn_country_codes_list = iucn_data{7};
    [threat_0, threat_1, threat_2] = cleanup_old_threats(iucn_data{4}, iucn_data{5}, iucn_data{6});
    processed_iucn_data.iucn_threats_list = combine_threats(threat_0, threat_1, threat_2, length(processed_iucn_data.q_row_list));

end


function [threat_0, threat_1, threat_2] = cleanup_old_threats(threat_0, threat_1, threat_2)
    
    threat_2(strcmp(threat_2, '1')) = {''};
    [threat_0, empties_0] = populate_empty_IDs(threat_0);
    [threat_1, empties_1] = populate_empty_IDs(threat_1);
    [threat_2, empties_2] = populate_empty_IDs(threat_2);

    threat_2_ind_test = cellfun(@(x) sum( char(x) == '.'), threat_2);
    threat_2_invalids = (threat_2_ind_test < 2) & (~empties_2);
    threat_1(threat_2_invalids) = threat_2(threat_2_invalids);

    threat_1_ind_test = cellfun(@(x) sum( char(x) == '.'), threat_1);
    threat_1_invalids = (threat_1_ind_test == 0) & (~empties_1);
    threat_0(threat_1_invalids) = threat_2(threat_1_invalids);

    threat_2(threat_2_invalids) = {'NaN'};
    threat_1(threat_1_invalids) = {'NaN'};
    

end


function [iucn_threats] = combine_threats(threat_0, threat_1, threat_2, iucn_data_length)
    iucn_threats = cell(iucn_data_length, 1);
    iucn_threats(~strcmp(threat_0, 'NaN')) = threat_0(~strcmp(threat_0, 'NaN'));
    iucn_threats(~strcmp(threat_1, 'NaN')) = threat_1(~strcmp(threat_1, 'NaN'));
    iucn_threats(~strcmp(threat_2, 'NaN')) = threat_2(~strcmp(threat_2, 'NaN'));
end


function threat_indexes = build_old_threat_indexes(iucn_data_length, threat_cause_class, threat_num, threat_0, threat_1, threat_2)

    threat_indexes = zeros(iucn_data_length, 1);
    for threat_ind = 1:threat_num %Populate threat vector by threat type e.g. populate many entries of threat type 2.3.1 simultaneously 
        current_threat = threat_cause_class(threat_ind);
        if (sum( char(current_threat) == '.') == 2)
            current_threat_set = strcmp(threat_2, current_threat);
        elseif (sum( char(current_threat) == '.') == 1)
            current_threat_set = strcmp(threat_1, current_threat);
        elseif (sum( char(current_threat) == '.') == 0)
            current_threat_set = strcmp(threat_0, current_threat);
        end

        threat_indexes(current_threat_set) = [threat_ind];
    end

end

function q_row = build_q_row(iucn_species_taxons)
    unique_species = unique(iucn_species_taxons);
    q_row = zeros(size(iucn_species_taxons));
        
    for species_ind = 1:length(unique_species)
        current_species_set = iucn_species_taxons == unique_species(species_ind);
        q_row(current_species_set) = species_ind;
    end
    
end


function [iucn_threats, iucn_lengths] = build_iucn_threat_data(iucn_threat_element, iucn_taxons, iucn_species_taxons, row_length)
    
    global_species_taxons = unique(iucn_species_taxons)';
    iucn_threats = cell(row_length, 1);
    
    for current_species = global_species_taxons
        current_taxon_indexes = (iucn_taxons == current_species);
        iucn_threats(iucn_species_taxons == current_species) = {iucn_threat_element(current_taxon_indexes)};
    end
    
    iucn_lengths = cellfun('length', iucn_threats);
    iucn_threats = vertcat(iucn_threats{:});
end

function expanded_iucn_data = expand_iucn_species_data(current_iucn_species_data_characteristic, row_length, iucn_lengths)
    expanded_iucn_data = cell(row_length, 1);
    for row_ind = 1:row_length
        expanded_iucn_data{row_ind} = repmat(current_iucn_species_data_characteristic(row_ind), [iucn_lengths(row_ind), 1]);
    end
    expanded_iucn_data = vertcat(expanded_iucn_data{:});
end


function [threat_cause_group] = generate_threat_cause_group(threat_cause_class)
    threat_cause_group = zeros(length(threat_cause_class), 1);
    threat_2_cond = cellfun(@(x) sum( char(x) == '.') == 2, threat_cause_class);
    threat_type_2_group = char(threat_cause_class(threat_2_cond));
    threat_cause_group(threat_2_cond) = str2num(threat_type_2_group(:, 1));
    threat_cause_group(~threat_2_cond) = floor(str2num(char(threat_cause_class(~threat_2_cond))));  
end


function ghg = simulate_ghg(ncoun, sector_num, scale_factor, build_type)
    ghg = cell(ncoun,1);
    for country_ind = 1:ncoun
        if strcmp(build_type, 'rand')
            ghg{country_ind} = scale_factor*rand(sector_num, 1);
        else ghg{country_ind} = scale_factor*ones(sector_num, 1);
        end
    end
end


function [eora_hscpc_concordance] = build_eora_hscpc_concordance(filepath, un_country_codes, iucn_country_codes, country_codes)
    
    eora_hscpc_concordance = cell(length(un_country_codes), 1);
    
    for country_ind = 1:length(un_country_codes)
        eora_hscpc_file = find_eora_to_hscpc(filepath, un_country_codes{country_ind}, 'i'); 
        if exist(eora_hscpc_file, 'file')
            eora2HS_raw = csvread(eora_hscpc_file);
            eora_hscpc_concordance{country_ind} = full(sparse(eora2HS_raw(:,2), eora2HS_raw(:,1),eora2HS_raw(:,3)));
        end
    end
    
    eora_hscpc_concordance = reorder_to_iucn(eora_hscpc_concordance, country_codes, iucn_country_codes);

end  


function iucn_object = reorder_to_iucn(object_to_reorder, country_codes, iucn_country_codes)

    iucn_country_code_length = length(country_codes); 
    [~, iucn_country_inds_to_use, un_country_inds_to_use] = intersect(country_codes, iucn_country_codes);
    
    iucn_object = cell(iucn_country_code_length, 1);
    iucn_object(iucn_country_inds_to_use) = object_to_reorder(un_country_inds_to_use);  
    
end


function threat_concordance = build_hscpc_threat_concordance(hscpc_concordance_filename, iucn_threat_characteristics, hscpc_sector_num)
    
    old_threat_cause_class_data = iucn_threat_characteristics.old_threat_cause_class;
    
    old_threat_cause_class = old_threat_cause_class_data(:, 1);
    concordance_data = csvread(hscpc_concordance_filename);
    old_threat_num = length(old_threat_cause_class_data(:, 1));
    old_threat_concordance = zeros(old_threat_num, hscpc_sector_num);
    
    for threat_class_ind = 1:old_threat_num;
        current_eco_sectors = concordance_data(concordance_data(:, 2) == threat_class_ind, 1);
        old_threat_concordance(threat_class_ind, current_eco_sectors) = true;
    end
    
    [~, new_threat_inds, old_threat_inds] = intersect(iucn_threat_characteristics.threat_cause_class, old_threat_cause_class); 
    threat_concordance = zeros(length(iucn_threat_characteristics.threat_num), hscpc_sector_num);
    threat_concordance(new_threat_inds, :) = old_threat_concordance(old_threat_inds, :);
    
end



function threat_concordance = build_eora_threat_concordance(eora_concordance_file_prefix, country_codes, un_to_iucn_codes, ncoun)
    
    ordered_un_codes = reorder_to_iucn(un_to_iucn_codes.un_country_codes, country_codes, un_to_iucn_codes.iucn_country_codes);
    threat_concordance = cell(ncoun, 1);
    for country_ind = 1:ncoun
		 concfile = [eora_concordance_file_prefix ordered_un_codes{country_ind} '=i.csv'];
         if ~exist(concfile)
            concfile = [eora_concordance_file_prefix ordered_un_codes{country_ind} '=c.csv'];
            if ~exist(concfile)
                concfile = [eora_concordance_file_prefix 'Ccc=i.csv'];
            end
         end
		 threat_concordance{country_ind} = dlmread(concfile); % Read concordance
    end
    
 end


function industry_indexes = build_hscpc_industry_indexes(iucn_data_length, threat_indexes_list, threat_concordance, threat_num)
    
    industry_indexes = cell(iucn_data_length, 1);    
    
    for threat_ind = 1:threat_num %Populate by threat type e.g. populate many entries of threat type 2.3.1 simultaneously
        current_C = find(threat_concordance(threat_ind, :));
        current_threat_set = (threat_indexes_list == threat_ind);
        industry_indexes(current_threat_set) = {current_C};
    end
    
end

function [industry_indexes] = build_eora_industry_indexes(iucn_data_length, threat_indexes_list, country_indexes_list, threat_concordance, ncoun)
    
    industry_indexes = cell(iucn_data_length, 1);    

    for country_ind = 1:ncoun
        current_country_indexes_list = (country_indexes_list == country_ind);               %select out current country
        current_threat_indexes_list = threat_indexes_list(current_country_indexes_list);    %select out threat indexes for that particular country
        current_threat_concordance = num2cell(threat_concordance{country_ind}, 2);          %select current threat concordance                 
        current_threat_concordance = cellfun(@(xx) find(xx), current_threat_concordance, 'UniformOutput', false);   %select indexes represented in concordance
        current_industries = current_threat_concordance(current_threat_indexes_list);       %
        industry_indexes(current_country_indexes_list) = current_industries;
    end

end


function processed_iucn_data = build_iucn_species_indexes(processed_iucn_data)
    
    disp('building index lists ...')

    processed_iucn_data.iucn_species_data.country_indexes_list = build_country_indexes_list(processed_iucn_data.iucn_species_data.iucn_data_length, processed_iucn_data.country_characteristics.ncoun, ...
                                                                            processed_iucn_data.country_characteristics.country_codes, processed_iucn_data.iucn_species_data.iucn_country_codes_list);
                                                                        
    processed_iucn_data.iucn_species_data.iucn_status_indexes_list = build_iucn_status_indexes_list(processed_iucn_data.iucn_species_data.iucn_data_length, length(processed_iucn_data.iucn_species_data.iucn_status_names), ...
                                                                                  processed_iucn_data.iucn_species_data.iucn_status_names, processed_iucn_data.iucn_species_data.iucn_status_list);
                                                                              
    processed_iucn_data.iucn_species_data.threat_indexes_list = build_threat_indexes_list(processed_iucn_data.iucn_species_data.iucn_data_length, processed_iucn_data.threat_characteristics.threat_cause_class, ...
                                                                        processed_iucn_data.threat_characteristics.threat_num, processed_iucn_data.iucn_species_data.iucn_threats_list);   

end


function country_indexes_list = build_country_indexes_list(iucn_data_length, ncoun, iucn_code_names, country_codes)
    
    country_indexes_list = zeros(iucn_data_length, 1);
    
    for country_ind = 1:ncoun
        current_country_ind_set = strcmp(iucn_code_names(country_ind), country_codes);
        country_indexes_list(current_country_ind_set) = country_ind;
    end

end

function iucn_status_indexes_list = build_iucn_status_indexes_list(iucn_data_length, iucn_status_num, threat_level_names, iucn_status)
    iucn_status_indexes_list = zeros(iucn_data_length, 1);

    for iucn_status_ind = 1:iucn_status_num %Populate threat vector by threat level type e.g. 'VU' etc

        current_threat_set = strcmp(iucn_status, threat_level_names(iucn_status_ind)); 
        iucn_status_indexes_list(current_threat_set) = iucn_status_ind;

    end

end


function threat_indexes_list = build_threat_indexes_list(iucn_data_length, threat_cause_class, threat_num, iucn_threat_list)
    threat_indexes_list = zeros(iucn_data_length, 1);
    for threat_ind = 1:threat_num %Populate threat vector by threat type e.g. populate many entries of threat type 2.3.1 simultaneously 
        current_threat = threat_cause_class(threat_ind);
        current_threat_set = strcmp(iucn_threat_list, current_threat);
        threat_indexes_list(current_threat_set) = [threat_ind];
    end

end



function scaled_industry_vals = scale_industry_vals(iucn_data_length, x, country_indexes_list, industry_indexes)
    scaled_industry_vals = cell(iucn_data_length, 1);

    for row_ind = 1:iucn_data_length
        
        current_x = x{country_indexes_list(row_ind)};
        
        if (~isempty(current_x))
            current_x = current_x(industry_indexes{row_ind});
            scaled_industry_vals{row_ind} = current_x./sum(current_x);
        end
        
    end

end


function [allCountriesCauseCounter, allCountriesRecordCounter] = build_allCountriesCauseCounter(nspec, NLEVELS, q_row_list, rl_indexes)
    
    allCountriesCauseCounter = zeros(nspec, NLEVELS);
    for row_ind = 1:length(q_row_list)
        sp = q_row_list(row_ind);
        status_level = rl_indexes(row_ind);
        allCountriesCauseCounter(sp,status_level) = allCountriesCauseCounter(sp,status_level) + 1;
    end

    allCountriesRecordCounter = allCountriesCauseCounter > 0;

end

