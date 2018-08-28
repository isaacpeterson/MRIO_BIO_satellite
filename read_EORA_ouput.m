addpath('/Users/E24661/Documents/MATLAB/BIO_SATELLITE/')
load_MRIO_objects = true;
country_of_interest = 'Australia';
status_levels_to_use = {'CR', 'EN', 'VU'}; %'all' or status levels
threat_direction = 'threats_from'; %'threats_on' or 'threats_from'
threat_type = 'international'; %'domestic' or 'international'
sort_type = 'species_num'; % 'threat_num' or 'species_num'
aggregate_type = 'IUCN_val'; % 'IUCN_val' or 'species_threat_path'
production_col = 2;
industry_col_to_use = 3;

if load_MRIO_objects == true
    datapath = '/Users/E24661/Documents/MATLAB/BIO_SATELLITE/EORA_outputs/';
    load('~/Documents/MATLAB/BIO_SATELLITE/IUCN_input_data/IUCN_data_object_for_manfred.mat')
    load([datapath 'SpThrSubs_domestic_final.mat'])
    load([datapath 'SpThrVals_domestic_final.mat'])
    load([datapath 'SpThrList_domestic_final.mat'])
    EORA_threat_tensor = sptensor(double(SpThrSubs), double(SpThrVals), double(max(SpThrSubs)));
    load([datapath 'satellite_species_characteristics.mat']);

    EORA_x_filename = [datapath, 'x_data_NCOUN_187.txt'];
    fid = fopen(EORA_x_filename);
        x_data = textscan(fid,'%s %s %s %s %f', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);
end

if (strcmp(status_levels_to_use, 'all'))
   species_to_use = 1:size(EORA_threat_tensor, 1);
else
    [~, ~, species_category_indexes] = intersect(species_characteristics.species_taxons, IUCN_data_object.IUCN_threat_taxons, 'stable');
    IUCN_threat_category = IUCN_data_object.IUCN_threat_status(species_category_indexes);
    species_to_use = find(ismember(IUCN_sort_category, status_levels_to_use)); 
end

industry_characteristics = struct();
industry_characteristics.country_names_list = x_data{1};
industry_characteristics.country_codes_list = x_data{2};
industry_characteristics.commodity_classification_list = x_data{4};

if strcmp(country_of_interest, 'all') 
    country_industries_to_use = 1:length(industry_characteristics.country_names_list);
    current_threat_tensor = EORA_threat_tensor(species_to_use, country_industries_to_use, :); 
else
    country_industries_to_use = find(strcmp(industry_characteristics.country_names_list, country_of_interest));
    if strcmp(threat_direction, 'threats_on')
        current_threat_tensor = EORA_threat_tensor(species_to_use, country_industries_to_use, :); 
    elseif strcmp(threat_direction, 'threats_from')
        current_threat_tensor = EORA_threat_tensor(species_to_use, :, country_industries_to_use); 
    end
end

[sorted_IUCN_vals, sorted_indexes] = sort(current_threat_tensor.vals, 'descend');
sorted_triplets = current_threat_tensor.subs(sorted_indexes, :);

%country specific species
sorted_species_taxons = sorted_triplets(:, 1); %directly read
sorted_species_names = species_characteristics.species_names(sorted_species_taxons);


sp_name = IUCN_data_object.IUCN_species_name;
[unique_sp, unique_sp_inds] = unique(IUCN_data_object.IUCN_species_name);
unique_sp_threat_num = zeros(numel(unique_sp), 1);

% build threat numbers from IUCN redlist
for i = 1:numel(unique_sp)
    unique_sp_threat_num(i) = sum(strcmp(sp_name, unique_sp(i)));
end

[sorted_IUCN_threat_num, sorted_indexes] = sort(unique_sp_threat_num, 'descend');
sorted_unique_IUCN_sp = unique_sp(sorted_indexes);

[unique_country_sp, ~, th_indexes] = intersect(sorted_species_names, sorted_unique_IUCN_sp, 'stable');

sorted_country_sp_threats = sorted_IUCN_threat_num(th_indexes);
threat_set = zeros(length(sorted_species_names), 1);

for current_ind = 1:length(th_indexes)
    current_set = strcmp(sorted_species_names, unique_country_sp(current_ind));
    threat_set(current_set) = sorted_IUCN_threat_num(th_indexes(current_ind));  
end












if strcmp(threat_direction, 'threats_on')
    
    sorted_production_industry_indexes = country_industries_to_use(sorted_triplets(:, production_col));  % subselection of industries means indexing goes from 1:N rather than the actual industry indexing
    sorted_consumption_industry_indexes = sorted_triplets(:, industry_col_to_use);
    
elseif strcmp(threat_direction, 'threats_from')
    sorted_production_industry_indexes = sorted_triplets(:, production_col);  % subselection of industries means indexing goes from 1:N rather than the actual industry indexing
    sorted_consumption_industry_indexes = country_industries_to_use(sorted_triplets(:, industry_col_to_use));    
end    


sorted_production_countries = industry_characteristics.country_names_list( sorted_production_industry_indexes );     
sorted_consumption_countries = industry_characteristics.country_names_list( sorted_consumption_industry_indexes );      
sorted_consumption_industries = industry_characteristics.commodity_classification_list(sorted_consumption_industry_indexes);
sorted_production_industries = industry_characteristics.commodity_classification_list(sorted_production_industry_indexes);
    

sorted_threat_cells = cell(1, 7);
sorted_threat_cells{1} = sorted_consumption_countries;
sorted_threat_cells{2} = sorted_consumption_industries;
sorted_threat_cells{3} = sorted_production_countries;
sorted_threat_cells{4} = sorted_production_industries;
sorted_threat_cells{5} = sorted_species_names;
sorted_threat_cells{6} = num2cell(threat_set);
sorted_threat_cells{7} = num2cell(sorted_IUCN_vals);

sorted_threat_array = [sorted_threat_cells{:}];

if strcmp(threat_direction, 'threats_on')
    domestic_threat_indexes = strcmp(sorted_consumption_countries, country_of_interest);
elseif strcmp(threat_direction, 'threats_from')
    domestic_threat_indexes = strcmp(sorted_production_countries, country_of_interest);
end
    
domestic_threat_array = sorted_threat_array(domestic_threat_indexes, :);
international_threat_array = sorted_threat_array(~domestic_threat_indexes, :);

if strcmp(threat_direction, 'threats_from')
    int_filename = strcat('~/Documents/MATLAB/BIO_SATELLITE/EORA_outputs/', threat_direction, '_', country_of_interest, '_on_', threat_type, '_ranked_by_threat_proportion.mat');
    dom_filename = strcat('~/Documents/MATLAB/BIO_SATELLITE/EORA_outputs/', threat_direction, '_', country_of_interest, '_on_', country_of_interest, '_ranked_by_threat_proportion.mat');
elseif  strcmp(threat_direction, 'threats_on')
    int_filename = strcat('~/Documents/MATLAB/BIO_SATELLITE/EORA_outputs/', threat_direction, '_', country_of_interest, '_from_', threat_type, '_ranked_by_threat_proportion.mat');
    dom_filename = strcat('~/Documents/MATLAB/BIO_SATELLITE/EORA_outputs/', threat_direction, '_', country_of_interest, '_from_', country_of_interest, '_ranked_by_threat_proportion.mat');
end

save(int_filename, 'international_threat_array', '-v7.3')
save(dom_filename, 'domestic_threat_array', '-v7.3')


industries_to_collate = sorted_consumption_industry_indexes;
unique_industries_to_collate = unique(industries_to_collate);

industries_to_aggregate = sorted_production_industry_indexes;
taxons_to_aggregate = sorted_species_taxons;
IUCN_vals_to_aggregate = sorted_IUCN_vals;

aggregated_species_counts = cell(length(unique_industries_to_collate), 1);

for current_ind = 1:length(unique_industries_to_collate)
    
    current_industry_ind = unique_industries_to_collate(current_ind);
    current_industry_set = (industries_to_collate == current_industry_ind);
    current_aggregate_set = industries_to_aggregate(current_industry_set);
    
    taxon_aggregate_set = taxons_to_aggregate(current_industry_set);
    IUCN_aggregate_set = IUCN_vals_to_aggregate(current_industry_set);
    
    [unique_aggregates, unique_inds] = unique(current_aggregate_set);
    
    if strcmp(aggregate_type, 'IUCN_val')
        current_IUCN_aggregates = zeros(length(unique_aggregates), 1);
        for agg_ind = 1:length(unique_aggregates)
            current_subset = current_aggregate_set == unique_aggregates(agg_ind);
            current_IUCN_aggregates(agg_ind) = sum(IUCN_aggregate_set(current_subset));
        end
    elseif strcmp(aggregate_type, 'species_threat_path')
        current_IUCN_aggregates = histc(current_aggregate_set, unique_aggregates);
    end
    aggregated_species_counts{current_industry_ind} = [repmat(current_industry_ind, size(unique_aggregates)), unique_aggregates, current_IUCN_aggregates];
    
end

aggregated_species_counts = vertcat(aggregated_species_counts{:});

[sorted_counts, sorted_inds] = sort(aggregated_species_counts(:, 3), 'descend');

sorted_species_counts = aggregated_species_counts(sorted_inds, :);
   
sorted_consumption_industry_indexes = sorted_species_counts(:, 1);
sorted_production_industry_indexes = sorted_species_counts(:, 2);
sorted_consumption_countries = industry_characteristics.country_names_list( sorted_consumption_industry_indexes);      
sorted_production_countries = industry_characteristics.country_names_list( sorted_production_industry_indexes );  
sorted_consumption_industries = industry_characteristics.commodity_classification_list(sorted_consumption_industry_indexes);
sorted_production_industries = industry_characteristics.commodity_classification_list(sorted_production_industry_indexes);
    
sorted_threat_cells = cell(1,6);
sorted_threat_cells{1} = sorted_consumption_countries;
sorted_threat_cells{2} = sorted_consumption_industries;
sorted_threat_cells{3} = sorted_production_countries;
sorted_threat_cells{4} = sorted_production_industries;
sorted_threat_cells{5} = num2cell(sorted_species_counts(:, 3));

%sorted_threat_cells{5} = species_characteristics.species_names(sorted_species_counts(:, 3));

%sorted_threat_cells{7} = num2cell(sorted_IUCN_vals);

sorted_threat_array = [sorted_threat_cells{:}];

save(strcat('~/Documents/MATLAB/BIO_SATELLITE/EORA_outputs/', threat_direction, '_', country_of_interest, '_ranked_by_industry_paths.mat'), 'sorted_threat_array', '-v7.3')




% if strcmp(threat_direction, 'threats_on')
%     domestic_threat_indexes = strcmp(sorted_consumption_countries, country_of_interest);
% elseif strcmp(threat_direction, 'threats_from')
%     domestic_threat_indexes = strcmp(sorted_production_countries, country_of_interest);
% end
%     
% domestic_threat_array = sorted_threat_array(domestic_threat_indexes, :);
% international_threat_array = sorted_threat_array(~domestic_threat_indexes, :);
% 
% if strcmp(threat_direction, 'threats_on')
%     industry_col_to_use = 4;
% elseif strcmp(threat_direction, 'threats_from')
%     industry_col_to_use = 2;
% end
% 
% if strcmp(threat_type, 'domestic')
% 	unique_product_industries = unique(domestic_threat_array(:, industry_col_to_use));
% elseif strcmp(threat_type, 'international')
% 	unique_product_industries = unique(international_threat_array(:, industry_col_to_use));  
% elseif strcmp(threat_type, 'all')
%     unique_product_industries = unique(sorted_threat_array(:, industry_col_to_use));
% end
%      
% industry_threat_rankings = zeros(length(unique_product_industries), 1);
% 
% for current_industry_ind = 1:length(unique_product_industries)
%     
%     if strcmp(threat_type, 'domestic')
%         current_industry_set = strcmp(domestic_threat_array(:, industry_col_to_use), unique_product_industries(current_industry_ind));
%     elseif strcmp(threat_type, 'international')
%         current_industry_set = strcmp(international_threat_array(:, industry_col_to_use), unique_product_industries(current_industry_ind));
%     elseif strcmp(threat_type, 'all')
%         current_industry_set = strcmp(sorted_threat_array(:, industry_col_to_use), unique_product_industries(current_industry_ind));
%     end
% 
%     if strcmp(sort_type, 'threat_num')
%         industry_threat_rankings(current_industry_ind) = sum(sorted_IUCN_vals(current_industry_set));
%     elseif strcmp(sort_type, 'species_num')
%         industry_threat_rankings(current_industry_ind) = length(unique(sorted_species_names(current_industry_set)));
%     end
%     
% end
% 
% [sorted_industry_ranks, sorted_indexes] = sort(industry_threat_rankings, 'descend');
% 
% industry_ranks = table(unique_product_industries(sorted_indexes), sorted_industry_ranks);
% 
% writetable(industry_ranks, strcat('~/Documents/MATLAB/BIO_SATELLITE/EORA_outputs/', country_of_interest, '_', threat_direction , '_', threat_type, '_ranked_by_', sort_type, '.xlsx'))


