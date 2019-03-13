function processed_iucn_data = process_raw_iucn_data(iucn_data_params)
    
    tic
    
    if (iucn_data_params.build_processed_iucn_data == false) && exist(iucn_data_params.processed_iucn_data_filename, 'file')
        disp(['loading processed iucn data from ', iucn_data_params.processed_iucn_data_filename])
        load([iucn_data_params.processed_iucn_data_filename])
    else
        disp(['processing iucn data to ', iucn_data_params.system_type, ' specification...'])
        if ~exist(iucn_data_params.processed_datapath, 'dir')
            mkdir(iucn_data_params.processed_datapath);
        end
        processed_iucn_data = build_processed_iucn_data(iucn_data_params);
    end

    disp([', processing and saving ', iucn_data_params.tensor_scale, ' level tensors...'])
    
    if (iucn_data_params.overwrite_tensors == true)
        build_iucn_tensors(processed_iucn_data, iucn_data_params);
    end
    
end

function processed_iucn_data = build_processed_iucn_data(iucn_data_params)
    
    un_to_iucn_codes = load_un_to_iucn_codes(iucn_data_params.un_to_iucn_codes_filename);
    eora_codes = load_eora_country_codes(iucn_data_params.eora_countries_filename);
    
    if strcmp(iucn_data_params.system_type, 'eora') 
        industry_characteristics = build_x_from_eora_data(iucn_data_params.eora_x_filename, un_to_iucn_codes); 
        [ghg_un, ghg_industry_codes] = build_eora_ghg(iucn_data_params.eora_ghg_filename, un_to_iucn_codes);    
    elseif strcmp(iucn_data_params.system_type, 'hscpc')
        industry_characteristics = build_x_from_hscpc_data(iucn_data_params.hscpc_x_filename, iucn_data_params.hscpc_country_codes_filename, un_to_iucn_codes);
    end
    
    
    if strcmp(iucn_data_params.iucn_data_type, 'new')  
       processed_iucn_data = expand_iucn_data(iucn_data_params.new_iucn_data_threats_filename, iucn_data_params.new_iucn_data_species_filename, ...
                                                    industry_characteristics.industry_codes_to_use, un_to_iucn_codes.iucn_country_codes);
    else
        processed_iucn_data = initialise_object_from_old_data(iucn_data_params.old_iucn_data_filename);
    end

    disp(['processed raw iucn data at ', num2str(toc)])
    
    processed_iucn_data = build_threat_classification_names(iucn_data_params.old_threat_cause_class_filename, iucn_data_params.read_threat_classification_from_file, processed_iucn_data);
    processed_iucn_data = build_iucn_country_names(iucn_data_params.read_iucn_countries_from_file, iucn_data_params.iucn_data_type, eora_codes, un_to_iucn_codes, processed_iucn_data);
        
    processed_iucn_data.un_to_iucn_codes = un_to_iucn_codes;
    processed_iucn_data.iucn_status_names = unique(processed_iucn_data.iucn_status_list); % 11 extinction levels (VUlnerable,EXtinct, etc)    
    processed_iucn_data.ncoun = length(processed_iucn_data.iucn_country_code_names);
    processed_iucn_data.iucn_status_num = length(processed_iucn_data.iucn_status_names);
    processed_iucn_data.nspec = max(processed_iucn_data.q_row_list);
    processed_iucn_data.iucn_data_length = length(processed_iucn_data.q_row_list);
    processed_iucn_data.threat_num = length(processed_iucn_data.threat_cause_class);
    
    disp(['building threat concordances ...'])
    if strcmp(iucn_data_params.system_type, 'hscpc')
        processed_iucn_data.threat_concordance = build_hscpc_threat_concordance(iucn_data_params.hscpc_concordance_filename, processed_iucn_data, iucn_data_params.hscpc_sector_num, processed_iucn_data.threat_cause_class, processed_iucn_data.threat_num);
    elseif strcmp(iucn_data_params.system_type, 'eora') 
        processed_iucn_data.threat_concordance = build_eora_threat_concordance(iucn_data_params.eora_concordance_file_prefix, processed_iucn_data.iucn_country_code_names, un_to_iucn_codes, processed_iucn_data.ncoun);
    end
   
    if strcmp(iucn_data_params.iucn_data_type, 'new')
        iucn_codes = un_to_iucn_codes.iucn_country_codes;
    else    
        iucn_codes = un_to_iucn_codes.un_country_codes;
    end

    processed_iucn_data.x = reorder_to_iucn(industry_characteristics.x_un, processed_iucn_data.iucn_country_code_names, iucn_codes);
    
    if strcmp(iucn_data_params.system_type, 'eora') 
        processed_iucn_data.x_names = reorder_to_iucn(industry_characteristics.x_un_names, processed_iucn_data.iucn_country_code_names, iucn_codes);
    elseif strcmp(iucn_data_params.system_type, 'hscpc')   
        
    end
    
    if (iucn_data_params.include_ghg == true)
        processed_iucn_data.ghg = reorder_to_iucn(ghg_un, processed_iucn_data.iucn_country_code_names, iucn_codes);
        processed_iucn_data.global_ghg = sum(cell2mat(processed_iucn_data.ghg)); 
    end
    
    processed_iucn_data.allcountriesflag = build_allcountriesflag(iucn_data_params.allcountriesflag_filename, processed_iucn_data.threat_cause_class, processed_iucn_data.old_threat_cause_class);
    processed_iucn_data = build_iucn_indexes(processed_iucn_data);
    
    processed_iucn_data.industry_indexes = build_industry_indexes(iucn_data_params.system_type, processed_iucn_data);
    processed_iucn_data.iucn_lengths = cellfun('length', processed_iucn_data.industry_indexes);
    
    processed_iucn_data.scaled_industry_vals = scale_industry_vals(processed_iucn_data.iucn_data_length, processed_iucn_data.x, processed_iucn_data.country_indexes_list, processed_iucn_data.industry_indexes);

    disp(['raw iucn data processed to ', iucn_data_params.system_type, ' specification at ' num2str(toc) ])
    
    if ~exist(iucn_data_params.processed_datapath, 'dir')
        mkdir(iucn_data_params.processed_datapath); 
    end
    
    if (iucn_data_params.save_processed_iucn_data == true)
        disp(['saving processed iucn data to ', iucn_data_params.processed_datapath, iucn_data_params.system_type, '/'])
        
        if ~exist([iucn_data_params.processed_datapath, iucn_data_params.system_type, '/'], 'dir')
            mkdir([iucn_data_params.processed_datapath, iucn_data_params.system_type, '/'])
        end
        
        save(iucn_data_params.processed_iucn_data_filename, 'processed_iucn_data', '-v7.3');
    end

    disp(['iucn data object built at ', toc])
    
end

function build_iucn_tensors(processed_iucn_data, iucn_data_params)

    disp(['writing tensors to ' iucn_data_params.pre_processed_tensor_filepath])
    
    if ~exist(iucn_data_params.pre_processed_tensor_filepath, 'dir')
        mkdir(iucn_data_params.pre_processed_tensor_filepath); 
    end
    
    if strcmp(iucn_data_params.tensor_scale, 'country')
        
        for country_index = 1:processed_iucn_data.ncoun
         
            rows_to_use = find(ismember(processed_iucn_data.country_indexes_list, country_index));
            current_iucn_tensor = build_current_tensor(processed_iucn_data, iucn_data_params.tensor_scale, rows_to_use);
            disp([processed_iucn_data.iucn_country_code_names{country_index}, ' tensor built at ' num2str(toc)])

            current_tensor_filename = [iucn_data_params.pre_processed_tensor_file_prefix, processed_iucn_data.iucn_country_code_names{country_index}, '.mat'];
            save(current_tensor_filename, 'current_iucn_tensor')
            
        end   
        
    else
        
        rows_to_use = 1:length(processed_iucn_data.country_indexes_list);
        current_iucn_tensor = build_current_tensor(processed_iucn_data, iucn_data_params.tensor_scale, rows_to_use);
        disp(['global iucn tensor built at ' num2str(toc)])
        	
        current_tensor_filename = [iucn_data_params.pre_processed_tensor_filepath, 'iucn_tensor.mat'];
        save(current_tensor_filename, 'current_iucn_tensor')
              
    end
    
end

function current_iucn_tensor = build_current_tensor(processed_iucn_data, tensor_type, rows_to_use)

    current_tensor_block = zeros(sum(processed_iucn_data.iucn_lengths(rows_to_use)), 1);
    if strcmp(tensor_type, 'global')
        current_tensor_block(:, 1) = build_iucn_tensor_indexes(processed_iucn_data.country_indexes_list(rows_to_use), processed_iucn_data.iucn_lengths(rows_to_use));
    else current_tensor_block(:, 1) = ones(sum(processed_iucn_data.iucn_lengths(rows_to_use)), 1);
    end
        
    current_tensor_block(:, 2) = build_iucn_tensor_indexes(processed_iucn_data.q_row_list(rows_to_use), processed_iucn_data.iucn_lengths(rows_to_use));
    current_tensor_block(:, 3) = cat(2, processed_iucn_data.industry_indexes{rows_to_use});
    current_tensor_block(:, 4) = build_iucn_tensor_indexes(processed_iucn_data.iucn_status_indexes_list(rows_to_use), processed_iucn_data.iucn_lengths(rows_to_use));
    current_tensor_block(:, 5) = build_iucn_tensor_indexes(processed_iucn_data.threat_indexes_list(rows_to_use), processed_iucn_data.iucn_lengths(rows_to_use));

    current_scaled_industry_vals = cat(1, processed_iucn_data.scaled_industry_vals{rows_to_use});
    current_iucn_tensor = sptensor(double(current_tensor_block), double(current_scaled_industry_vals), double(max(current_tensor_block)));
   
end

function tensor_index_array = build_iucn_tensor_indexes(current_processed_iucn_data, current_lengths)

    tensor_index_array = arrayfun(@(x, y) repmat(x, [y, 1]), current_processed_iucn_data, current_lengths, 'UniformOutput', false);
    tensor_index_array = cell2mat(tensor_index_array); 
    
end

function allcountriesflag_new = build_allcountriesflag(allcountriesflag_filename, new_threat_cause_class, old_threat_cause_class)
    load(allcountriesflag_filename)
    [~, old_inds, new_inds] = intersect(old_threat_cause_class, new_threat_cause_class);
    allcountriesflag_new = zeros(size(new_threat_cause_class));
    allcountriesflag_new(new_inds) = allcountriesflag(old_inds);  
end


function processed_iucn_data = build_iucn_country_names(read_iucn_countries_from_file, iucn_data_type, eora_codes, un_to_iucn_codes, processed_iucn_data)
    
    if read_iucn_countries_from_file
       
       if strcmp(iucn_data_type, 'old') 
            
            iucn_country_names = eora_codes.eora_country_names;
            iucn_country_code_names = eora_codes.eora_country_codes;
            empty_names = strcmp(iucn_country_code_names, ' ');
           
       elseif strcmp(iucn_data_type, 'new')
           
            [iucn_country_code_names, sort_inds] = sort(un_to_iucn_codes.iucn_country_codes);
            iucn_country_names = un_to_iucn_codes.iucn_country_names(sort_inds);
            empty_names = cellfun('isempty', iucn_country_code_names);
       
       end
       
   else 
        [iucn_country_code_names, unique_inds, ~] = unique(processed_iucn_data.iucn_country_codes_list);
        iucn_country_names = processed_iucn_data.country_names_list(unique_inds); 
        empty_names = cellfun('isempty', iucn_country_code_names);

   end

   processed_iucn_data.iucn_country_code_names = iucn_country_code_names(~empty_names);
   processed_iucn_data.iucn_country_names = iucn_country_names(~empty_names);
   
end


function processed_iucn_data = build_threat_classification_names(old_threat_cause_class_filename, read_threat_classification_from_file, processed_iucn_data)
    
    disp(['building threat classifications ...'])
    
    fid = fopen(old_threat_cause_class_filename);
        old_threat_cause_class_data = textscan(fid,'%s %s %f', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);

    processed_iucn_data.old_threat_cause_class = strrep(old_threat_cause_class_data{1}, '''','');
    processed_iucn_data.old_threat_cause_names = old_threat_cause_class_data{2};
    processed_iucn_data.old_threat_cause_group = old_threat_cause_class_data{3};
        
    if read_threat_classification_from_file
        processed_iucn_data.threat_cause_class = processed_iucn_data.old_threat_cause_class;
        processed_iucn_data.threat_cause_group = processed_iucn_data.old_threat_cause_group;
   else
       processed_iucn_data.threat_cause_class = unique(processed_iucn_data.iucn_threats_list);

       [~, name_inds, old_name_inds] = intersect(processed_iucn_data.threat_cause_class, processed_iucn_data.old_threat_cause_class, 'stable');
       threat_cause_names = cell(length(processed_iucn_data.threat_cause_class), 1);
       threat_cause_names(name_inds) = processed_iucn_data.old_threat_cause_names(old_name_inds);
       unlisted_set = setdiff(1:length(processed_iucn_data.threat_cause_class), name_inds);
       threat_cause_names(unlisted_set) = processed_iucn_data.threat_cause_class(unlisted_set);
       
       processed_iucn_data.threat_cause_names = threat_cause_names;
       processed_iucn_data.threat_cause_group = generate_threat_cause_group(processed_iucn_data.threat_cause_class);
    end 
   
end

%un_to_iucn_codes_filename = iucn_data_params.un_to_iucn_codes_filename;
function [un_to_iucn_codes] = load_un_to_iucn_codes(un_to_iucn_codes_filename)
    fid = fopen(un_to_iucn_codes_filename);
        un_to_iucn_data = textscan(fid,'%s %s %s', 'delimiter', ';');
    fclose(fid);
    un_to_iucn_codes = struct();
    un_to_iucn_codes.iucn_country_codes = un_to_iucn_data{3};
    un_to_iucn_codes.un_country_codes = un_to_iucn_data{2};
    un_to_iucn_codes.iucn_country_names = un_to_iucn_data{1};
    
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



function [industry_characteristics] = build_x_from_eora_data(eora_x_filename, un_to_iucn_codes)       

    fid = fopen(eora_x_filename);
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
    
    [industry_characteristics.x_un, industry_characteristics.industry_codes_to_use] = build_iucn_industries_data(industry_characteristics.numerical_industry_data, industry_characteristics.country_codes_list, industry_data_list, commodities_data_list, un_to_iucn_codes);
                                    
    [industry_characteristics.x_un_names, ~] = build_iucn_industries_data(industry_characteristics.commodity_classification_list, industry_characteristics.country_codes_list, industry_data_list, commodities_data_list, un_to_iucn_codes); 
                                    
end  

function [industry_characteristics] = build_x_from_hscpc_data(hscpc_x_filename, hscpc_country_codes_filename, un_to_iucn_codes)
    
    un_country_codes = un_to_iucn_codes.un_country_codes;
    fid = fopen(hscpc_country_codes_filename);
        hscpc_country_list = textscan(fid, '%s %s', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);
    
    load(hscpc_x_filename)
    [~, ~, un_country_inds] = intersect(hscpc_country_list{2}, un_country_codes); 
    
    industry_characteristics = struct();
    country_num = length(un_country_codes);
    industry_characteristics.x_un = cell(country_num, 1);
     
    for country_ind = 1:country_num
        industry_characteristics.x_un{un_country_inds(country_ind)} = GlobalRoot(:, country_ind);   
    end
    
    industry_characteristics.industry_codes_to_use = ~cellfun('isempty', industry_characteristics.x_un);
end

% function hscpc_x = build_x_from_hscpc_data_b(filepath, un_to_iucn_codes, iucn_country_code_names, ncoun)
% 
%     load([filepath 'Global_Root.mat'])
%     hscpc_x = cell(ncoun, 1);
%     
%     [~, ~, country_inds_to_use] = intersect(iucn_country_code_names, un_to_iucn_codes{3}); 
%     
%     for country_ind = 1:length(country_inds_to_use);
%         hscpc_x{country_ind} = GlobalRoot(:, country_inds_to_use(country_ind));   
%     end
%     
% end

% new_iucn_data_threats_filename = iucn_data_params.new_iucn_data_threats_filename;
% new_iucn_data_species_filename = iucn_data_params.new_iucn_data_species_filename;
% industry_codes_to_use;
% iucn_country_codes = un_to_iucn_codes.iucn_country_codes;
                                                
function processed_iucn_data = expand_iucn_data(new_iucn_data_threats_filename, new_iucn_data_species_filename, industry_codes_to_use, iucn_country_codes)

    fid = fopen(new_iucn_data_threats_filename);
    iucn_threat_data = textscan(fid,'%f %f %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);

    fid = fopen(new_iucn_data_species_filename);
        iucn_species_data = textscan(fid,'%f %f %s %s %s %s %s %s %s %s %s %s %s %s %s', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);

    iucn_threats = iucn_threat_data{14};
    empty_threats = cellfun('isempty', iucn_threats); 
    
    iucn_threats = iucn_threats(~empty_threats); %remove empty threats from threat data
    iucn_threat_taxons = iucn_threat_data{2}(~empty_threats);
    iucn_threat_status = iucn_threat_data{10}(~empty_threats);

    iucn_data_country_codes = iucn_species_data{12};
    iucn_data_country_names = iucn_species_data{13};
    iucn_data_species_taxons = iucn_species_data{2};
    
    [iucn_data_country_codes, sorted_iucn_code_inds] = sort(iucn_data_country_codes);
    iucn_data_country_names = iucn_data_country_names(sorted_iucn_code_inds);
    iucn_data_species_taxons = iucn_data_species_taxons(sorted_iucn_code_inds);
    
    iucn_codes_to_use = unique(iucn_country_codes(industry_codes_to_use));  
    
    iucn_codes_to_use = iucn_codes_to_use(~cellfun('isempty', iucn_codes_to_use));  %remove empty code names from UN list consequently removes empties from iucn list
                                                                                    
    iucn_inds_to_use = ismember(iucn_data_country_codes, iucn_codes_to_use); %find all intersecting countries in iucn_list and UN list, i.e. remove iucn countries that there is no industrial data for.

    iucn_data_country_codes = iucn_data_country_codes(iucn_inds_to_use);
    iucn_data_species_taxons = iucn_data_species_taxons(iucn_inds_to_use);

    excluded_countries = unique(iucn_data_country_names(~iucn_inds_to_use));
    iucn_data_country_names = iucn_data_country_names(iucn_inds_to_use);

    row_length = length(iucn_data_country_codes);

    [expanded_iucn_threats, iucn_lengths] = build_iucn_threat_data(iucn_threats, iucn_threat_taxons, iucn_data_species_taxons, row_length);    
    [expanded_iucn_status, ~] = build_iucn_threat_data(iucn_threat_status, iucn_threat_taxons, iucn_data_species_taxons, row_length);
    expanded_iucn_status = strrep(expanded_iucn_status, '/', '_');
    expanded_iucn_data_taxons = expand_iucn_species_data(row_length, iucn_data_species_taxons, iucn_lengths);
    expanded_country_codes = expand_iucn_species_data(row_length, iucn_data_country_codes, iucn_lengths);
    expanded_country_names = expand_iucn_species_data(row_length, iucn_data_country_names, iucn_lengths);
    expanded_q_row = build_q_row(expanded_iucn_data_taxons);
    
    processed_iucn_data = struct;
    processed_iucn_data.excluded_countries = excluded_countries;
    processed_iucn_data.country_names_list = expanded_country_names;
    processed_iucn_data.iucn_threats_list = expanded_iucn_threats;
    processed_iucn_data.iucn_status_list = expanded_iucn_status;
    processed_iucn_data.iucn_taxons_list = expanded_iucn_data_taxons;
    processed_iucn_data.iucn_country_codes_list = expanded_country_codes;
    processed_iucn_data.q_row_list = expanded_q_row;
    
    processed_iucn_data.iucn_threat_taxons = iucn_threat_taxons;
    processed_iucn_data.iucn_threat_status = iucn_threat_status;
    
    processed_iucn_data.iucn_species_kingdom = iucn_threat_data{3}(~empty_threats);
    processed_iucn_data.iucn_species_phylum = iucn_threat_data{4}(~empty_threats);
    processed_iucn_data.iucn_species_class = iucn_threat_data{5}(~empty_threats);
    processed_iucn_data.iucn_species_order = iucn_threat_data{6}(~empty_threats);
    processed_iucn_data.iucn_species_family = iucn_threat_data{7}(~empty_threats);
    processed_iucn_data.iucn_species_genus = iucn_threat_data{8}(~empty_threats);
    processed_iucn_data.iucn_species_name = iucn_threat_data{9}(~empty_threats);
    
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

%build_iucn_threat_data(iucn_threats_by_taxon, iucn_taxons, iucn_data_species_taxons, row_length);   

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

function expanded_iucn_data = expand_iucn_species_data(row_length, iucn_species_element, iucn_lengths)
    expanded_iucn_data = cell(row_length, 1);
    for row_ind = 1:row_length
        expanded_iucn_data{row_ind} = repmat(iucn_species_element(row_ind), [iucn_lengths(row_ind), 1]);
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


 
function [ghg_un, ghg_industry_codes] = build_eora_ghg(eora_ghg_filename, un_to_iucn_codes)
    
    fid = fopen(eora_ghg_filename);
        ghg_eora_Data = textscan(fid,'%f %s %s %s %s %f', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);
    
    industry_characteristics = struct();
    industry_characteristics.country_codes_list = ghg_eora_Data{3};
    industry_data_list = strcmp(ghg_eora_Data{4}, 'Industries');
    commodities_data_list = strcmp(ghg_eora_Data{4}, 'Commodities');
    industry_characteristics.numerical_industry_data = ghg_eora_Data{6};
    
    industry_names = ghg_eora_Data{5};
    
%     export_data = regexpi(industry_names,'export');
%     non_export_data = cellfun('isempty', export_data);
    
    [ghg_un, ghg_industry_codes] = build_iucn_industries_data(industry_characteristics.numerical_industry_data, industry_characteristics.country_codes_list, ...
                                                                industry_data_list, commodities_data_list, un_to_iucn_codes);
                                    
%     export_data = regexpi(industry_names,'export');
%     non_export_data = cellfun('isempty', export_data);
%     
%     ghg_industries = Build_eora_ghg(industry_data, ghg_country_codes_list, un_country_codes, ghg_CO2);
%     ghg_comms = Build_eora_ghg(commodities_data, ghg_country_codes_list, un_country_codes, ghg_CO2);
%     
%     empty_industries = cellfun('isempty', ghg_industries);
%     empty_comms = cellfun('isempty', ghg_comms);
%     
%     inds_to_replace = empty_industries & ~empty_comms;
%     %data_to_use = industry_data & (ghg_CO2 > 0);
%     eora_ghg = ghg_industries;
%     eora_ghg(inds_to_replace) = ghg_comms(inds_to_replace);
    
end

                                                                                                            
function [industries_data_un, industry_codes_to_use] = build_iucn_industries_data(data_to_sort, country_codes_list, industry_data_list, commodities_data_list, un_to_iucn_codes)

    industries_data_un = sort_eora_data_to_un(country_codes_list(industry_data_list), data_to_sort(industry_data_list), un_to_iucn_codes.un_country_codes);
    commodities_data_un = sort_eora_data_to_un(country_codes_list(commodities_data_list), data_to_sort(industry_data_list), un_to_iucn_codes.un_country_codes);
    
    empty_industries = cellfun('isempty', industries_data_un);
    empty_comms = cellfun('isempty', commodities_data_un);
    
    inds_to_replace = empty_industries & ~empty_comms;

    industries_data_un(inds_to_replace) = commodities_data_un(inds_to_replace);
    industry_codes_to_use = ~cellfun('isempty', industries_data_un);
    
end

function sorted_eora_data = sort_eora_data_to_un(country_codes_list_to_use, data_to_use, un_country_codes)

    data_code_names = unique(country_codes_list_to_use);

    [~, ~, un_country_inds] = intersect(data_code_names, un_country_codes);
    
    sorted_eora_data = cell(length(un_country_codes), 1);
    
    for country_ind = 1:length(data_code_names)
        current_set = strcmp(country_codes_list_to_use, data_code_names(country_ind));
        sorted_eora_data{un_country_inds(country_ind)} = data_to_use(current_set);   
    end
    
end

% function eora_ghg = Build_eora_ghg(data_to_use, ghg_country_codes_list, un_country_codes, ghg_CO2)
% 
%     ghg_CO2_to_use = ghg_CO2(data_to_use);
%     ghg_country_codes_list_to_use = ghg_country_codes_list(data_to_use);
%     eora_code_names = unique(ghg_country_codes_list_to_use);
% 
%     [~, ~, un_country_inds] = intersect(eora_code_names, un_country_codes);
%     
%     eora_ghg = cell(length(un_country_codes), 1);
%     
%     for country_ind = 1:length(eora_code_names)
%         current_set = strcmp(ghg_country_codes_list_to_use, eora_code_names(country_ind));
%         eora_ghg{un_country_inds(country_ind)} = ghg_CO2_to_use(current_set);    %trim off re-export re-import for each one
%     end
%     
% end
%filepath = '~/Documents/MATLAB/BIO_SATELLITE/hscpc_to_eora_concs/'
%un_country_codes = un_to_iucn_codes{2};
%build_eora_hscpc_Concordance('~/Documents/MATLAB/BIO_SATELLITE/hscpc_to_eora_concs/', processed_iucn_data.un_to_iucn_codes{2}, processed_iucn_data.un_to_iucn_codes{3}, processed_iucn_data.iucn_country_code_names)

function [eora_hscpc_concordance] = build_eora_hscpc_concordance(filepath, un_country_codes, iucn_country_codes, iucn_country_code_names)
    
    eora_hscpc_concordance = cell(length(un_country_codes), 1);
    
    for country_ind = 1:length(un_country_codes)
        eora_hscpc_file = find_eora_to_hscpc(filepath, un_country_codes{country_ind}, 'i'); 
        if exist(eora_hscpc_file, 'file')
            eora2HS_raw = csvread(eora_hscpc_file);
            eora_hscpc_concordance{country_ind} = full(sparse(eora2HS_raw(:,2), eora2HS_raw(:,1),eora2HS_raw(:,3)));
        end
    end
    
    eora_hscpc_concordance = reorder_to_iucn(eora_hscpc_concordance, iucn_country_code_names, iucn_country_codes);

end  


% object_to_reorder = x;
% iucn_country_codes = un_to_iucn_codes.iucn_country_codes;

function iucn_object = reorder_to_iucn(object_to_reorder, iucn_country_code_names, iucn_country_codes)

    iucn_country_code_length = length(iucn_country_code_names); 
    [~, iucn_country_inds_to_use, un_country_inds_to_use] = intersect(iucn_country_code_names, iucn_country_codes);
    
    iucn_object = cell(iucn_country_code_length, 1);
    iucn_object(iucn_country_inds_to_use) = object_to_reorder(un_country_inds_to_use);  
    
end



function threat_concordance = build_hscpc_threat_concordance(hscpc_concordance_filename, processed_iucn_data, hscpc_sector_num, threat_cause_class, threat_num)
    
    old_threat_cause_class_data = processed_iucn_data.old_threat_cause_class;
    
    old_threat_cause_class = old_threat_cause_class_data(:, 1);
    concordance_data = csvread(hscpc_concordance_filename);
    old_threat_num = length(old_threat_cause_class_data(:, 1));
    old_threat_concordance = zeros(old_threat_num, hscpc_sector_num);
    
    for threat_class_ind = 1:old_threat_num;
        current_eco_sectors = concordance_data(concordance_data(:, 2) == threat_class_ind, 1);
        old_threat_concordance(threat_class_ind, current_eco_sectors) = true;
    end
    
    [~, new_threat_inds, old_threat_inds] = intersect(threat_cause_class, old_threat_cause_class); 
    threat_concordance = zeros(length(threat_num), hscpc_sector_num);
    threat_concordance(new_threat_inds, :) = old_threat_concordance(old_threat_inds, :);
    
end



function threat_concordance = build_eora_threat_concordance(eora_concordance_file_prefix, iucn_country_code_names, un_to_iucn_codes, ncoun)
    
    ordered_un_codes = reorder_to_iucn(un_to_iucn_codes.un_country_codes, iucn_country_code_names, un_to_iucn_codes.iucn_country_codes);
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


function [old_threat_cause_class, old_threat_class_data] = load_old_threat_cause_class(filenm) 
    old_threat_cause_class = cellread(filenm);
    old_threat_class_data =  strrep(old_threat_cause_class,'''',''); %Remove the ' characters used in the CSV file for Excel-safety
end

function [industry_indexes] = build_industry_indexes(system_type, processed_iucn_data)
   
   if strcmp(system_type, 'eora')
       industry_indexes = build_eora_industry_indexes(processed_iucn_data.iucn_data_length, processed_iucn_data.threat_indexes_list, processed_iucn_data.country_indexes_list, processed_iucn_data.threat_concordance, processed_iucn_data.ncoun);
   elseif strcmp(system_type, 'hscpc')
       industry_indexes = build_hscpc_industry_indexes(processed_iucn_data.iucn_data_length, processed_iucn_data.threat_indexes_list, processed_iucn_data.threat_concordance, processed_iucn_data.threat_num);
   end
   
end

function [industry_indexes] = build_hscpc_industry_indexes(iucn_data_length, threat_indexes_list, threat_concordance, threat_num)
    
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


function processed_iucn_data = build_iucn_indexes(processed_iucn_data)
    disp(['building index lists ...'])
    processed_iucn_data.country_indexes_list = build_country_indexes_list(processed_iucn_data.iucn_data_length, processed_iucn_data.ncoun, processed_iucn_data.iucn_country_code_names, processed_iucn_data.iucn_country_codes_list);
    processed_iucn_data.iucn_status_indexes_list = build_iucn_status_indexes_list(processed_iucn_data.iucn_data_length, processed_iucn_data.iucn_status_num, processed_iucn_data.iucn_status_names, processed_iucn_data.iucn_status_list);
    processed_iucn_data.threat_indexes_list = build_threat_indexes_list(processed_iucn_data.iucn_data_length, processed_iucn_data.threat_cause_class, processed_iucn_data.threat_num, processed_iucn_data.iucn_threats_list);   
    processed_iucn_data.data_length = length(processed_iucn_data.country_indexes_list);
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

%%%% CC routines finished

function [allCountriesCauseCounter, allCountriesRecordCounter] = build_allCountriesCauseCounter(nspec, NLEVELS, q_row_list, rl_indexes)
    
    allCountriesCauseCounter = zeros(nspec, NLEVELS);
    for row_ind = 1:length(q_row_list)
        sp = q_row_list(row_ind);
        status_level = rl_indexes(row_ind);
        allCountriesCauseCounter(sp,status_level) = allCountriesCauseCounter(sp,status_level) + 1;
    end

    allCountriesRecordCounter = allCountriesCauseCounter > 0;

end

