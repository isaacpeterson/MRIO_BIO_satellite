function IUCN_data_object = process_IUCN_data_routines(satellite_params)

    if (satellite_params.build_IUCN_data_object == false) && exist([satellite_params.output_data_filepath, satellite_params.IUCN_data_object_filename], 'file')
        disp(['loading processed IUCN data from ', satellite_params.output_data_filepath, satellite_params.IUCN_data_object_filename])
        load([satellite_params.output_data_filepath, satellite_params.IUCN_data_object_filename])
        return
    else
        disp(['processing IUCN data to ', satellite_params.system_type, ' specification...'])
        if ~exist(satellite_params.output_data_filepath, 'dir')
            mkdir(satellite_params.output_data_filepath);
        end
    end

    UN_to_IUCN_codes = load_UN_to_IUCN_codes(satellite_params.UN_to_IUCN_codes_filename);
    EORA_codes = load_EORA_country_codes(satellite_params.EORA_countries_filename);
    
    if strcmp(satellite_params.system_type, 'EORA') 
        industry_characteristics = build_x_from_EORA_data(satellite_params.EORA_x_filename, UN_to_IUCN_codes); 
        [GHG_UN, GHG_industry_codes] = build_EORA_GHG(satellite_params.EORA_GHG_filename, UN_to_IUCN_codes);    
    elseif strcmp(satellite_params.system_type, 'HSCPC')
        industry_characteristics = build_x_from_HSCPC_data(satellite_params.HSCPC_x_filename, satellite_params.HSCPC_country_codes_filename, UN_to_IUCN_codes);
    end
    
    
    
    if strcmp(satellite_params.IUCN_data_type, 'new')  
       IUCN_data_object = process_raw_IUCN_data(satellite_params.new_IUCN_data_threats_filename, satellite_params.new_IUCN_data_species_filename, ...
                                                    industry_characteristics.industry_codes_to_use, UN_to_IUCN_codes.IUCN_industry_codes);
    else
        IUCN_data_object = initialise_object_from_old_data(satellite_params.old_IUCN_data_filename);
    end

    disp(['processed raw IUCN data at ', num2str(toc)])
    
    IUCN_data_object = build_threat_classification_names(satellite_params.old_threat_cause_class_filename, satellite_params.read_threat_classification_from_file, IUCN_data_object);
    IUCN_data_object = build_IUCN_country_names(satellite_params.read_IUCN_countries_from_file, satellite_params.IUCN_data_type, EORA_codes, UN_to_IUCN_codes, IUCN_data_object);
        
    IUCN_data_object.UN_to_IUCN_codes = UN_to_IUCN_codes;
    IUCN_data_object.IUCN_status_names = unique(IUCN_data_object.IUCN_status_list); % 11 extinction levels (VUlnerable,EXtinct, etc)    
    IUCN_data_object.NCOUN = length(IUCN_data_object.IUCN_country_code_names);
    IUCN_data_object.IUCN_status_num = length(IUCN_data_object.IUCN_status_names);
    IUCN_data_object.NSPEC = max(IUCN_data_object.q_row_list);
    IUCN_data_object.IUCN_data_length = length(IUCN_data_object.q_row_list);
    IUCN_data_object.threat_num = length(IUCN_data_object.threat_cause_class);
    
    disp(['building threat concordances ...'])
    if strcmp(satellite_params.system_type, 'HSCPC')
        IUCN_data_object.threat_concordance = build_HSCPC_threat_concordance(satellite_params.HSCPC_concordance_filename, IUCN_data_object, satellite_params.HSCPC_sector_num, IUCN_data_object.threat_cause_class, IUCN_data_object.threat_num);
    elseif strcmp(satellite_params.system_type, 'EORA') 
        IUCN_data_object.threat_concordance = build_EORA_threat_concordance(satellite_params.EORA_concordance_file_prefix, IUCN_data_object.IUCN_country_code_names, UN_to_IUCN_codes, IUCN_data_object.NCOUN);
    end
   
    if strcmp(satellite_params.IUCN_data_type, 'new')
        IUCN_codes = UN_to_IUCN_codes.IUCN_industry_codes;
    else    
        IUCN_codes = UN_to_IUCN_codes.UN_industry_codes;
    end

    IUCN_data_object.x = reorder_to_IUCN(industry_characteristics.x_UN, IUCN_data_object.IUCN_country_code_names, IUCN_codes);
    
    if strcmp(satellite_params.system_type, 'EORA') 
        IUCN_data_object.x_names = reorder_to_IUCN(industry_characteristics.x_UN_names, IUCN_data_object.IUCN_country_code_names, IUCN_codes);
    elseif strcmp(satellite_params.system_type, 'HSCPC')   
        
    end
    
    if (satellite_params.include_GHG == true)
        IUCN_data_object.GHG = reorder_to_IUCN(GHG_UN, IUCN_data_object.IUCN_country_code_names, IUCN_codes);
        IUCN_data_object.global_GHG = sum(cell2mat(IUCN_data_object.GHG)); 
    end
    
    
    IUCN_data_object.allcountriesflag = build_allcountriesflag(satellite_params.allcountriesflag_filename, IUCN_data_object.threat_cause_class, IUCN_data_object.old_threat_cause_class);
    IUCN_data_object = build_IUCN_indexes(IUCN_data_object);
    
    IUCN_data_object.industry_indexes = build_industry_indexes(satellite_params.system_type, IUCN_data_object);
    IUCN_data_object.IUCN_lengths = cellfun('length', IUCN_data_object.industry_indexes);
    
    IUCN_data_object.scaled_industry_vals = scale_industry_vals(IUCN_data_object.IUCN_data_length, IUCN_data_object.x, IUCN_data_object.country_indexes_list, IUCN_data_object.industry_indexes);

    disp(['IUCN data object processed to ', satellite_params.system_type, ' specification at ' num2str(toc) ])
    
    if ~exist(satellite_params.output_data_filepath, 'dir')
        mkdir(satellite_params.output_data_filepath); 
    end
    
    if (satellite_params.save_processed_IUCN_data == true)
        disp(['saving processed IUCN data to ', satellite_params.IUCN_data_object_filename])
        save([satellite_params.output_data_filepath, satellite_params.IUCN_data_object_filename], 'IUCN_data_object', '-v7.3');
    end

    disp(['IUCN data object built at ' toc ', processing and saving ', satellite_params.tensor_scale, ' level tensors...'])
    
    build_IUCN_tensors(IUCN_data_object, satellite_params);
    
end


function build_IUCN_tensors(IUCN_data_object, IUCN_data_params)

    disp(['writing tensors to ' IUCN_data_params.tensor_folder])
    
    if ~exist(IUCN_data_params.tensor_folder, 'dir')
        mkdir(IUCN_data_params.tensor_folder); 
    end
    
    if strcmp(IUCN_data_params.tensor_scale, 'country')
        
        for country_index = 1:IUCN_data_object.NCOUN
         
            rows_to_use = find(ismember(IUCN_data_object.country_indexes_list, country_index));
            current_IUCN_tensor = build_current_tensor(IUCN_data_object, IUCN_data_params.tensor_scale, rows_to_use);
            disp([IUCN_data_object.IUCN_country_code_names{country_index}, ' tensor built at ' num2str(toc)])
            if (IUCN_data_params.save_IUCN_tensors == true)
            	current_tensor_filename = [IUCN_data_params.tensor_folder, 'IUCN_tensor_', IUCN_data_object.IUCN_country_code_names{country_index}, '.mat'];
            	save(current_tensor_filename, 'current_IUCN_tensor')
            end
            
        end   
        
    else
        
        rows_to_use = 1:length(IUCN_data_object.country_indexes_list);
        current_IUCN_tensor = build_current_tensor(IUCN_data_object, IUCN_data_params.tensor_scale, rows_to_use);
        disp(['global IUCN tensor built at ' num2str(toc)])
        
        if (IUCN_data_params.save_IUCN_tensors == true)
        	current_tensor_filename = [IUCN_data_params.tensor_folder, 'IUCN_tensor.mat'];
            save(current_tensor_filename, 'current_IUCN_tensor')
        end
              
    end
    
end

function current_IUCN_tensor = build_current_tensor(IUCN_data_object, tensor_type, rows_to_use)

    current_tensor_block = zeros(sum(IUCN_data_object.IUCN_lengths(rows_to_use)), 1);
    if strcmp(tensor_type, 'global')
        current_tensor_block(:, 1) = build_IUCN_tensor_indexes(IUCN_data_object.country_indexes_list(rows_to_use), IUCN_data_object.IUCN_lengths(rows_to_use));
    else current_tensor_block(:, 1) = ones(sum(IUCN_data_object.IUCN_lengths(rows_to_use)), 1);
    end
        
    current_tensor_block(:, 2) = build_IUCN_tensor_indexes(IUCN_data_object.q_row_list(rows_to_use), IUCN_data_object.IUCN_lengths(rows_to_use));
    current_tensor_block(:, 3) = cat(2, IUCN_data_object.industry_indexes{rows_to_use});
    current_tensor_block(:, 4) = build_IUCN_tensor_indexes(IUCN_data_object.IUCN_status_indexes_list(rows_to_use), IUCN_data_object.IUCN_lengths(rows_to_use));
    current_tensor_block(:, 5) = build_IUCN_tensor_indexes(IUCN_data_object.threat_indexes_list(rows_to_use), IUCN_data_object.IUCN_lengths(rows_to_use));

    current_scaled_industry_vals = cat(1, IUCN_data_object.scaled_industry_vals{rows_to_use});
    current_IUCN_tensor = sptensor(double(current_tensor_block), double(current_scaled_industry_vals), double(max(current_tensor_block)));
   
end

function tensor_index_array = build_IUCN_tensor_indexes(current_IUCN_data_object, current_lengths)

    tensor_index_array = arrayfun(@(x, y) repmat(x, [y, 1]), current_IUCN_data_object, current_lengths, 'UniformOutput', false);
    tensor_index_array = cell2mat(tensor_index_array); 
    
end

function allcountriesflag_new = build_allcountriesflag(allcountriesflag_filename, new_threat_cause_class, old_threat_cause_class)
    load(allcountriesflag_filename)
    [~, old_inds, new_inds] = intersect(old_threat_cause_class, new_threat_cause_class);
    allcountriesflag_new = zeros(size(new_threat_cause_class));
    allcountriesflag_new(new_inds) = allcountriesflag(old_inds);  
end


function IUCN_data_object = build_IUCN_country_names(read_IUCN_countries_from_file, IUCN_data_type, EORA_codes, UN_to_IUCN_codes, IUCN_data_object)
    
    if read_IUCN_countries_from_file
       
       if strcmp(IUCN_data_type, 'old') 
            
            IUCN_country_names = EORA_codes.EORA_country_names;
            IUCN_country_code_names = EORA_codes.EORA_country_codes;
            empty_names = strcmp(IUCN_country_code_names, ' ');
           
       elseif strcmp(IUCN_data_type, 'new')
           
            [IUCN_country_code_names, sort_inds] = sort(UN_to_IUCN_codes.IUCN_industry_codes);
            IUCN_country_names = UN_to_IUCN_codes.IUCN_country_names(sort_inds);
            empty_names = cellfun('isempty', IUCN_country_code_names);
       
       end
       
   else 
        [IUCN_country_code_names, unique_inds, ~] = unique(IUCN_data_object.IUCN_country_codes_list);
        IUCN_country_names = IUCN_data_object.country_names_list(unique_inds); 
        empty_names = cellfun('isempty', IUCN_country_code_names);

   end

   IUCN_data_object.IUCN_country_code_names = IUCN_country_code_names(~empty_names);
   IUCN_data_object.IUCN_country_names = IUCN_country_names(~empty_names);
   
end


function IUCN_data_object = build_threat_classification_names(old_threat_cause_class_filename, read_threat_classification_from_file, IUCN_data_object)
    
    disp(['building threat classifications ...'])
    
    fid = fopen(old_threat_cause_class_filename);
        old_threat_cause_class_data = textscan(fid,'%s %s %f', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);

    IUCN_data_object.old_threat_cause_class = strrep(old_threat_cause_class_data{1}, '''','');
    IUCN_data_object.old_threat_cause_names = old_threat_cause_class_data{2};
    IUCN_data_object.old_threat_cause_group = old_threat_cause_class_data{3};
        
    if read_threat_classification_from_file
        IUCN_data_object.threat_cause_class = IUCN_data_object.old_threat_cause_class;
        IUCN_data_object.threat_cause_group = IUCN_data_object.old_threat_cause_group;
   else
       IUCN_data_object.threat_cause_class = unique(IUCN_data_object.IUCN_threats_list);

       [~, name_inds, old_name_inds] = intersect(IUCN_data_object.threat_cause_class, IUCN_data_object.old_threat_cause_class, 'stable');
       threat_cause_names = cell(length(IUCN_data_object.threat_cause_class), 1);
       threat_cause_names(name_inds) = IUCN_data_object.old_threat_cause_names(old_name_inds);
       unlisted_set = setdiff(1:length(IUCN_data_object.threat_cause_class), name_inds);
       threat_cause_names(unlisted_set) = IUCN_data_object.threat_cause_class(unlisted_set);
       
       IUCN_data_object.threat_cause_names = threat_cause_names;
       IUCN_data_object.threat_cause_group = generate_threat_cause_group(IUCN_data_object.threat_cause_class);
    end 
   
end

%UN_to_IUCN_codes_filename = IUCN_data_params.UN_to_IUCN_codes_filename;
function [UN_to_IUCN_codes] = load_UN_to_IUCN_codes(UN_to_IUCN_codes_filename)
    fid = fopen(UN_to_IUCN_codes_filename);
        UN_to_IUCN_data = textscan(fid,'%s %s %s', 'delimiter', ';');
    fclose(fid);
    UN_to_IUCN_codes = struct();
    UN_to_IUCN_codes.IUCN_industry_codes = UN_to_IUCN_data{3};
    UN_to_IUCN_codes.UN_industry_codes = UN_to_IUCN_data{2};
    UN_to_IUCN_codes.IUCN_country_names = UN_to_IUCN_data{1};
    
end
    
function [EORA_codes] = load_EORA_country_codes(EORA_countries_filename)
    [~,~, EORA_country_data] = xlsread(EORA_countries_filename); 
    EORA_country_names = EORA_country_data(3:end,3);
    EORA_country_codes = EORA_country_data(3:end,4);
    EORA_country_codes = strrep(EORA_country_codes, 'TZM', 'TZA');   %BIG NOTE - THIS WAS THE OTHER WAY AROUND IN THE PREVIOUS CODE - i.e. TANZANIA WAS EXCLUDED FROM THE ANALYSIS
    empty_names = strcmp(EORA_country_codes, ' ');
    EORA_country_codes = EORA_country_codes(~empty_names);
    EORA_country_names = EORA_country_names(~empty_names);
    [EORA_codes.EORA_country_codes, sorted_inds] = sort(EORA_country_codes);
    EORA_codes.EORA_country_names = EORA_country_names(sorted_inds);
end   



function [industry_characteristics] = build_x_from_EORA_data(EORA_x_filename, UN_to_IUCN_codes)       

    fid = fopen(EORA_x_filename);
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
    industry_characteristics.numerical_industry_data = x_data{5};
        
    industry_data_list = strcmp(x_data{3}, 'Industries');
    commodities_data_list = strcmp(x_data{3}, 'Commodities');
    
    [industry_characteristics.x_UN, industry_characteristics.industry_codes_to_use] = build_IUCN_industries_data(industry_characteristics.numerical_industry_data, industry_characteristics.country_codes_list, industry_data_list, commodities_data_list, UN_to_IUCN_codes);
                                    
    [industry_characteristics.x_UN_names, ~] = build_IUCN_industries_data(industry_characteristics.commodity_classification_list, industry_characteristics.country_codes_list, industry_data_list, commodities_data_list, UN_to_IUCN_codes); 
                                    
end  

function [industry_characteristics] = build_x_from_HSCPC_data(HSCPC_x_filename, HSCPC_country_codes_filename, UN_to_IUCN_codes)
    
    UN_industry_codes = UN_to_IUCN_codes.UN_industry_codes;
    fid = fopen(HSCPC_country_codes_filename);
        HSCPC_country_list = textscan(fid, '%s %s', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);
    
    load(HSCPC_x_filename)
    [~, ~, UN_country_inds] = intersect(HSCPC_country_list{2}, UN_industry_codes); 
    
    industry_characteristics = struct();
    country_num = length(UN_industry_codes);
    industry_characteristics.x_UN = cell(country_num, 1);
     
    for country_ind = 1:country_num
        industry_characteristics.x_UN{UN_country_inds(country_ind)} = GlobalRoot(:, country_ind);   
    end
    
    industry_characteristics.industry_codes_to_use = ~cellfun('isempty', industry_characteristics.x_UN);
end

% function HSCPC_x = build_x_from_HSCPC_data_b(filepath, UN_to_IUCN_codes, IUCN_country_code_names, NCOUN)
% 
%     load([filepath 'Global_Root.mat'])
%     HSCPC_x = cell(NCOUN, 1);
%     
%     [~, ~, country_inds_to_use] = intersect(IUCN_country_code_names, UN_to_IUCN_codes{3}); 
%     
%     for country_ind = 1:length(country_inds_to_use);
%         HSCPC_x{country_ind} = GlobalRoot(:, country_inds_to_use(country_ind));   
%     end
%     
% end

% new_IUCN_data_threats_filename = IUCN_data_params.new_IUCN_data_threats_filename;
% new_IUCN_data_species_filename = IUCN_data_params.new_IUCN_data_species_filename;
% industry_codes_to_use;
% IUCN_industry_codes = UN_to_IUCN_codes.IUCN_industry_codes;
                                                
function IUCN_data_object = process_raw_IUCN_data(new_IUCN_data_threats_filename, new_IUCN_data_species_filename, industry_codes_to_use, IUCN_industry_codes)

    fid = fopen(new_IUCN_data_threats_filename);
    IUCN_threat_data = textscan(fid,'%f %f %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);

    fid = fopen(new_IUCN_data_species_filename);
        IUCN_species_data = textscan(fid,'%f %f %s %s %s %s %s %s %s %s %s %s %s %s %s', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);

    IUCN_threats = IUCN_threat_data{14};
    empty_threats = cellfun('isempty', IUCN_threats); 
    
    IUCN_threats = IUCN_threats(~empty_threats); %remove empty threats from threat data
    IUCN_threat_taxons = IUCN_threat_data{2}(~empty_threats);
    IUCN_threat_status = IUCN_threat_data{10}(~empty_threats);

    IUCN_data_country_codes = IUCN_species_data{12};
    IUCN_data_country_names = IUCN_species_data{13};
    IUCN_data_species_taxons = IUCN_species_data{2};
    
    [IUCN_data_country_codes, sorted_IUCN_code_inds] = sort(IUCN_data_country_codes);
    IUCN_data_country_names = IUCN_data_country_names(sorted_IUCN_code_inds);
    IUCN_data_species_taxons = IUCN_data_species_taxons(sorted_IUCN_code_inds);
    
    IUCN_codes_to_use = unique(IUCN_industry_codes(industry_codes_to_use));  
    
    IUCN_codes_to_use = IUCN_codes_to_use(~cellfun('isempty', IUCN_codes_to_use));  %remove empty code names from UN list consequently removes empties from IUCN list
                                                                                    
    IUCN_inds_to_use = ismember(IUCN_data_country_codes, IUCN_codes_to_use); %find all intersecting countries in IUCN_list and UN list, i.e. remove IUCN countries that there is no industrial data for.

    IUCN_data_country_codes = IUCN_data_country_codes(IUCN_inds_to_use);
    IUCN_data_species_taxons = IUCN_data_species_taxons(IUCN_inds_to_use);

    excluded_countries = unique(IUCN_data_country_names(~IUCN_inds_to_use));
    IUCN_data_country_names = IUCN_data_country_names(IUCN_inds_to_use);

    row_length = length(IUCN_data_country_codes);

    [expanded_IUCN_threats, IUCN_lengths] = build_IUCN_threat_data(IUCN_threats, IUCN_threat_taxons, IUCN_data_species_taxons, row_length);    
    [expanded_IUCN_status, ~] = build_IUCN_threat_data(IUCN_threat_status, IUCN_threat_taxons, IUCN_data_species_taxons, row_length);
    expanded_IUCN_status = strrep(expanded_IUCN_status, '/', '_');
    expanded_IUCN_data_taxons = expand_IUCN_species_data(row_length, IUCN_data_species_taxons, IUCN_lengths);
    expanded_country_codes = expand_IUCN_species_data(row_length, IUCN_data_country_codes, IUCN_lengths);
    expanded_country_names = expand_IUCN_species_data(row_length, IUCN_data_country_names, IUCN_lengths);
    expanded_q_row = build_q_row(expanded_IUCN_data_taxons);
    
    IUCN_data_object = struct;
    IUCN_data_object.excluded_countries = excluded_countries;
    IUCN_data_object.country_names_list = expanded_country_names;
    IUCN_data_object.IUCN_threats_list = expanded_IUCN_threats;
    IUCN_data_object.IUCN_status_list = expanded_IUCN_status;
    IUCN_data_object.IUCN_taxons_list = expanded_IUCN_data_taxons;
    IUCN_data_object.IUCN_country_codes_list = expanded_country_codes;
    IUCN_data_object.q_row_list = expanded_q_row;
    
    IUCN_data_object.IUCN_threat_taxons = IUCN_threat_taxons;
    IUCN_data_object.IUCN_threat_status = IUCN_threat_status;
    
    IUCN_data_object.IUCN_species_kingdom = IUCN_threat_data{3}(~empty_threats);
    IUCN_data_object.IUCN_species_phylum = IUCN_threat_data{4}(~empty_threats);
    IUCN_data_object.IUCN_species_class = IUCN_threat_data{5}(~empty_threats);
    IUCN_data_object.IUCN_species_order = IUCN_threat_data{6}(~empty_threats);
    IUCN_data_object.IUCN_species_family = IUCN_threat_data{7}(~empty_threats);
    IUCN_data_object.IUCN_species_genus = IUCN_threat_data{8}(~empty_threats);
    IUCN_data_object.IUCN_species_name = IUCN_threat_data{9}(~empty_threats);
    
end


function IUCN_data_object = initialise_object_from_old_data(old_IUCN_data_filename)
    fid = fopen(old_IUCN_data_filename);
        IUCN_data = textscan(fid,'%f %s %s %s %s %s %s', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);

    IUCN_data_object.q_row_list = IUCN_data{1};
    IUCN_data_object.country_names_list = IUCN_data{2};
    IUCN_data_object.IUCN_status_list = IUCN_data{3}; % 11 threat levels (VUlnerable,EXtinct, etc)
    IUCN_data_object.IUCN_country_codes_list = IUCN_data{7};
    [threat_0, threat_1, threat_2] = cleanup_old_threats(IUCN_data{4}, IUCN_data{5}, IUCN_data{6});
    IUCN_data_object.IUCN_threats_list = combine_threats(threat_0, threat_1, threat_2, length(IUCN_data_object.q_row_list));

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


function [IUCN_threats] = combine_threats(threat_0, threat_1, threat_2, IUCN_data_length)
    IUCN_threats = cell(IUCN_data_length, 1);
    IUCN_threats(~strcmp(threat_0, 'NaN')) = threat_0(~strcmp(threat_0, 'NaN'));
    IUCN_threats(~strcmp(threat_1, 'NaN')) = threat_1(~strcmp(threat_1, 'NaN'));
    IUCN_threats(~strcmp(threat_2, 'NaN')) = threat_2(~strcmp(threat_2, 'NaN'));
end


function threat_indexes = build_old_threat_indexes(IUCN_data_length, threat_cause_class, threat_num, threat_0, threat_1, threat_2)

    threat_indexes = zeros(IUCN_data_length, 1);
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

function q_row = build_q_row(IUCN_species_taxons)
    unique_species = unique(IUCN_species_taxons);
    q_row = zeros(size(IUCN_species_taxons));
        
    for species_ind = 1:length(unique_species)
        current_species_set = IUCN_species_taxons == unique_species(species_ind);
        q_row(current_species_set) = species_ind;
    end
    
end

%build_IUCN_threat_data(IUCN_threats_by_taxon, IUCN_taxons, IUCN_data_species_taxons, row_length);   

function [IUCN_threats, IUCN_lengths] = build_IUCN_threat_data(IUCN_threat_element, IUCN_taxons, IUCN_species_taxons, row_length)
    
    global_species_taxons = unique(IUCN_species_taxons)';
    IUCN_threats = cell(row_length, 1);
    for current_species = global_species_taxons
        current_taxon_indexes = (IUCN_taxons == current_species);
        IUCN_threats(IUCN_species_taxons == current_species) = {IUCN_threat_element(current_taxon_indexes)};
    end
    IUCN_lengths = cellfun('length', IUCN_threats);
    IUCN_threats = vertcat(IUCN_threats{:});
end

function expanded_IUCN_data = expand_IUCN_species_data(row_length, IUCN_species_element, IUCN_lengths)
    expanded_IUCN_data = cell(row_length, 1);
    for row_ind = 1:row_length
        expanded_IUCN_data{row_ind} = repmat(IUCN_species_element(row_ind), [IUCN_lengths(row_ind), 1]);
    end
    expanded_IUCN_data = vertcat(expanded_IUCN_data{:});
end


function [threat_cause_group] = generate_threat_cause_group(threat_cause_class)
    threat_cause_group = zeros(length(threat_cause_class), 1);
    threat_2_cond = cellfun(@(x) sum( char(x) == '.') == 2, threat_cause_class);
    threat_type_2_group = char(threat_cause_class(threat_2_cond));
    threat_cause_group(threat_2_cond) = str2num(threat_type_2_group(:, 1));
    threat_cause_group(~threat_2_cond) = floor(str2num(char(threat_cause_class(~threat_2_cond))));  
end




function GHG = simulate_GHG(NCOUN, sector_num, scale_factor, build_type)
    GHG = cell(NCOUN,1);
    for country_ind = 1:NCOUN
        if strcmp(build_type, 'rand')
            GHG{country_ind} = scale_factor*rand(sector_num, 1);
        else GHG{country_ind} = scale_factor*ones(sector_num, 1);
        end
    end
end


 
function [GHG_UN, GHG_industry_codes] = build_EORA_GHG(EORA_GHG_filename, UN_to_IUCN_codes)
    
    fid = fopen(EORA_GHG_filename);
        GHG_EORA_Data = textscan(fid,'%f %s %s %s %s %f', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);
    
    industry_characteristics = struct();
    industry_characteristics.country_codes_list = GHG_EORA_Data{3};
    industry_data_list = strcmp(GHG_EORA_Data{4}, 'Industries');
    commodities_data_list = strcmp(GHG_EORA_Data{4}, 'Commodities');
    industry_characteristics.numerical_industry_data = GHG_EORA_Data{6};
    
    industry_names = GHG_EORA_Data{5};
    
%     export_data = regexpi(industry_names,'export');
%     non_export_data = cellfun('isempty', export_data);
    
    [GHG_UN, GHG_industry_codes] = build_IUCN_industries_data(industry_characteristics.numerical_industry_data, industry_characteristics.country_codes_list, ...
                                                                industry_data_list, commodities_data_list, UN_to_IUCN_codes);
                                    
%     export_data = regexpi(industry_names,'export');
%     non_export_data = cellfun('isempty', export_data);
%     
%     GHG_industries = Build_EORA_GHG(industry_data, GHG_country_codes_list, UN_industry_codes, GHG_CO2);
%     GHG_comms = Build_EORA_GHG(commodities_data, GHG_country_codes_list, UN_industry_codes, GHG_CO2);
%     
%     empty_industries = cellfun('isempty', GHG_industries);
%     empty_comms = cellfun('isempty', GHG_comms);
%     
%     inds_to_replace = empty_industries & ~empty_comms;
%     %data_to_use = industry_data & (GHG_CO2 > 0);
%     EORA_GHG = GHG_industries;
%     EORA_GHG(inds_to_replace) = GHG_comms(inds_to_replace);
    
end

                                                                                                            
function [industries_data_UN, industry_codes_to_use] = build_IUCN_industries_data(data_to_sort, country_codes_list, industry_data_list, commodities_data_list, UN_to_IUCN_codes)

    industries_data_UN = sort_EORA_data_to_UN(country_codes_list(industry_data_list), data_to_sort(industry_data_list), UN_to_IUCN_codes.UN_industry_codes);
    commodities_data_UN = sort_EORA_data_to_UN(country_codes_list(commodities_data_list), data_to_sort(industry_data_list), UN_to_IUCN_codes.UN_industry_codes);
    
    empty_industries = cellfun('isempty', industries_data_UN);
    empty_comms = cellfun('isempty', commodities_data_UN);
    
    inds_to_replace = empty_industries & ~empty_comms;

    industries_data_UN(inds_to_replace) = commodities_data_UN(inds_to_replace);
    industry_codes_to_use = ~cellfun('isempty', industries_data_UN);
    
end

function sorted_EORA_data = sort_EORA_data_to_UN(country_codes_list_to_use, data_to_use, UN_industry_codes)

    data_code_names = unique(country_codes_list_to_use);

    [~, ~, UN_country_inds] = intersect(data_code_names, UN_industry_codes);
    
    sorted_EORA_data = cell(length(UN_industry_codes), 1);
    
    for country_ind = 1:length(data_code_names)
        current_set = strcmp(country_codes_list_to_use, data_code_names(country_ind));
        sorted_EORA_data{UN_country_inds(country_ind)} = data_to_use(current_set);   
    end
    
end

% function EORA_GHG = Build_EORA_GHG(data_to_use, GHG_country_codes_list, UN_industry_codes, GHG_CO2)
% 
%     GHG_CO2_to_use = GHG_CO2(data_to_use);
%     GHG_country_codes_list_to_use = GHG_country_codes_list(data_to_use);
%     EORA_code_names = unique(GHG_country_codes_list_to_use);
% 
%     [~, ~, UN_country_inds] = intersect(EORA_code_names, UN_industry_codes);
%     
%     EORA_GHG = cell(length(UN_industry_codes), 1);
%     
%     for country_ind = 1:length(EORA_code_names)
%         current_set = strcmp(GHG_country_codes_list_to_use, EORA_code_names(country_ind));
%         EORA_GHG{UN_country_inds(country_ind)} = GHG_CO2_to_use(current_set);    %trim off re-export re-import for each one
%     end
%     
% end
%filepath = '~/Documents/MATLAB/BIO_SATELLITE/HSCPC_to_EORA_concs/'
%UN_industry_codes = UN_to_IUCN_codes{2};
%build_EORA_HSCPC_Concordance('~/Documents/MATLAB/BIO_SATELLITE/HSCPC_to_EORA_concs/', IUCN_data_object.UN_to_IUCN_codes{2}, IUCN_data_object.UN_to_IUCN_codes{3}, IUCN_data_object.IUCN_country_code_names)

function [EORA_HSCPC_Concordance] = build_EORA_HSCPC_Concordance(filepath, UN_industry_codes, IUCN_industry_codes, IUCN_country_code_names)
    
    EORA_HSCPC_concordance = cell(length(UN_industry_codes), 1);
    
    for country_ind = 1:length(UN_industry_codes)
        EORA_HSCPC_file = find_EORA_to_HSCPC(filepath, UN_industry_codes{country_ind}, 'i'); 
        if exist(EORA_HSCPC_file)
            EORA2HS_raw = csvread(EORA_HSCPC_file);
            EORA_HSCPC_concordance{country_ind} = full(sparse(EORA2HS_raw(:,2), EORA2HS_raw(:,1),EORA2HS_raw(:,3)));
        end
    end
    
    EORA_HSCPC_concordance = reorder_to_IUCN(EORA_HSCPC_concordance, IUCN_country_code_names, IUCN_industry_codes);

end  


% object_to_reorder = x;
% IUCN_industry_codes = UN_to_IUCN_codes.IUCN_industry_codes;

function IUCN_object = reorder_to_IUCN(object_to_reorder, IUCN_country_code_names, IUCN_industry_codes)

    IUCN_country_code_length = length(IUCN_country_code_names); 
    [~, IUCN_country_inds_to_use, UN_country_inds_to_use] = intersect(IUCN_country_code_names, IUCN_industry_codes);
    
    IUCN_object = cell(IUCN_country_code_length, 1);
    IUCN_object(IUCN_country_inds_to_use) = object_to_reorder(UN_country_inds_to_use);  
    
end



function threat_concordance = build_HSCPC_threat_concordance(HSCPC_concordance_filename, IUCN_data_object, HSCPC_sector_num, threat_cause_class, threat_num)
    
    old_threat_cause_class_data = IUCN_data_object.old_threat_cause_class;
    
    old_threat_cause_class = old_threat_cause_class_data(:, 1);
    concordance_data = csvread(HSCPC_concordance_filename);
    old_threat_num = length(old_threat_cause_class_data(:, 1));
    old_threat_concordance = zeros(old_threat_num, HSCPC_sector_num);
    
    for threat_class_ind = 1:old_threat_num;
        current_eco_sectors = concordance_data(concordance_data(:, 2) == threat_class_ind, 1);
        old_threat_concordance(threat_class_ind, current_eco_sectors) = true;
    end
    
    [~, new_threat_inds, old_threat_inds] = intersect(threat_cause_class, old_threat_cause_class); 
    threat_concordance = zeros(length(threat_num), HSCPC_sector_num);
    threat_concordance(new_threat_inds, :) = old_threat_concordance(old_threat_inds, :);
    
end



function threat_concordance = build_EORA_threat_concordance(EORA_concordance_file_prefix, IUCN_country_code_names, UN_to_IUCN_codes, NCOUN)
    
    ordered_UN_codes = reorder_to_IUCN(UN_to_IUCN_codes.UN_industry_codes, IUCN_country_code_names, UN_to_IUCN_codes.IUCN_industry_codes);
    threat_concordance = cell(NCOUN, 1);
    for country_ind = 1:NCOUN
		 ConcFile = [EORA_concordance_file_prefix ordered_UN_codes{country_ind} '=i.csv'];
         if ~exist(ConcFile)
            ConcFile = [EORA_concordance_file_prefix ordered_UN_codes{country_ind} '=c.csv'];
            if ~exist(ConcFile)
                ConcFile = [EORA_concordance_file_prefix 'Ccc=i.csv'];
            end
         end
		 threat_concordance{country_ind} = dlmread(ConcFile); % Read concordance
    end
    
 end


function [old_threat_cause_class, old_threat_class_data] = load_old_threat_cause_class(filenm) 
    old_threat_cause_class = cellread(filenm);
    old_threat_class_data =  strrep(old_threat_cause_class,'''',''); %Remove the ' characters used in the CSV file for Excel-safety
end

function [industry_indexes] = build_industry_indexes(system_type, IUCN_data_object)
   
   if strcmp(system_type, 'EORA')
       industry_indexes = build_EORA_industry_indexes(IUCN_data_object.IUCN_data_length, IUCN_data_object.threat_indexes_list, IUCN_data_object.country_indexes_list, IUCN_data_object.threat_concordance, IUCN_data_object.NCOUN);
   elseif strcmp(system_type, 'HSCPC')
       industry_indexes = build_HSCPC_industry_indexes(IUCN_data_object.IUCN_data_length, IUCN_data_object.threat_indexes_list, IUCN_data_object.threat_concordance, IUCN_data_object.threat_num);
   end
   
end

function [industry_indexes] = build_HSCPC_industry_indexes(IUCN_data_length, threat_indexes_list, threat_concordance, threat_num)
    
    industry_indexes = cell(IUCN_data_length, 1);    
    
    for threat_ind = 1:threat_num %Populate by threat type e.g. populate many entries of threat type 2.3.1 simultaneously
        current_C = find(threat_concordance(threat_ind, :));
        current_threat_set = (threat_indexes_list == threat_ind);
        industry_indexes(current_threat_set) = {current_C};
    end
    
end

function [industry_indexes] = build_EORA_industry_indexes(IUCN_data_length, threat_indexes_list, country_indexes_list, threat_concordance, NCOUN)
    
    industry_indexes = cell(IUCN_data_length, 1);    

    for country_ind = 1:NCOUN
        current_country_indexes_list = (country_indexes_list == country_ind);               %select out current country
        current_threat_indexes_list = threat_indexes_list(current_country_indexes_list);    %select out threat indexes for that particular country
        current_threat_concordance = num2cell(threat_concordance{country_ind}, 2);          %select current threat concordance                 
        current_threat_concordance = cellfun(@(xx) find(xx), current_threat_concordance, 'UniformOutput', false);   %select indexes represented in concordance
        current_industries = current_threat_concordance(current_threat_indexes_list);       %
        industry_indexes(current_country_indexes_list) = current_industries;
    end

end


function IUCN_data_object = build_IUCN_indexes(IUCN_data_object)
    disp(['building index lists ...'])
    IUCN_data_object.country_indexes_list = build_country_indexes_list(IUCN_data_object.IUCN_data_length, IUCN_data_object.NCOUN, IUCN_data_object.IUCN_country_code_names, IUCN_data_object.IUCN_country_codes_list);
    IUCN_data_object.IUCN_status_indexes_list = build_IUCN_status_indexes_list(IUCN_data_object.IUCN_data_length, IUCN_data_object.IUCN_status_num, IUCN_data_object.IUCN_status_names, IUCN_data_object.IUCN_status_list);
    IUCN_data_object.threat_indexes_list = build_threat_indexes_list(IUCN_data_object.IUCN_data_length, IUCN_data_object.threat_cause_class, IUCN_data_object.threat_num, IUCN_data_object.IUCN_threats_list);   
    IUCN_data_object.data_length = length(IUCN_data_object.country_indexes_list);
end


function country_indexes_list = build_country_indexes_list(IUCN_data_length, NCOUN, IUCN_code_names, country_codes)
    country_indexes_list = zeros(IUCN_data_length, 1);
    
    for country_ind = 1:NCOUN
        current_country_ind_set = strcmp(IUCN_code_names(country_ind), country_codes);
        country_indexes_list(current_country_ind_set) = country_ind;
    end

end

function IUCN_status_indexes_list = build_IUCN_status_indexes_list(IUCN_data_length, IUCN_status_num, threat_level_names, IUCN_status)
    IUCN_status_indexes_list = zeros(IUCN_data_length, 1);

    for IUCN_status_ind = 1:IUCN_status_num %Populate threat vector by threat level type e.g. 'VU' etc

        current_threat_set = strcmp(IUCN_status, threat_level_names(IUCN_status_ind)); 
        IUCN_status_indexes_list(current_threat_set) = IUCN_status_ind;

    end

end


function threat_indexes_list = build_threat_indexes_list(IUCN_data_length, threat_cause_class, threat_num, IUCN_threat_list)
    threat_indexes_list = zeros(IUCN_data_length, 1);
    for threat_ind = 1:threat_num %Populate threat vector by threat type e.g. populate many entries of threat type 2.3.1 simultaneously 
        current_threat = threat_cause_class(threat_ind);
        current_threat_set = strcmp(IUCN_threat_list, current_threat);
        threat_indexes_list(current_threat_set) = [threat_ind];
    end

end



function scaled_industry_vals = scale_industry_vals(IUCN_data_length, x, country_indexes_list, industry_indexes)
    scaled_industry_vals = cell(IUCN_data_length, 1);

    for row_ind = 1:IUCN_data_length
        
        current_x = x{country_indexes_list(row_ind)};
        
        if (~isempty(current_x))
            current_x = current_x(industry_indexes{row_ind});
            scaled_industry_vals{row_ind} = current_x./sum(current_x);
        end
        
    end

end

%%%% CC routines finished

function [allCountriesCauseCounter, allCountriesRecordCounter] = build_allCountriesCauseCounter(NSPEC, NLEVELS, q_row_list, rl_indexes)
    
    allCountriesCauseCounter = zeros(NSPEC, NLEVELS);
    for row_ind = 1:length(q_row_list)
        sp = q_row_list(row_ind);
        status_level = rl_indexes(row_ind);
        allCountriesCauseCounter(sp,status_level) = allCountriesCauseCounter(sp,status_level) + 1;
    end

    allCountriesRecordCounter = allCountriesCauseCounter > 0;

end
