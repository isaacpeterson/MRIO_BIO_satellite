function footprint_objects = build_input_objects(input_type, footprint_objects_params, system_type, satellite_species_characteristics)
    
    if strcmp(input_type, 'industry')
    
        footprint_objects.un_to_iucn_codes = load_un_to_iucn_codes(footprint_objects_params.un_to_iucn_codes_filename);
    
        if strcmp(system_type, 'eora') 
            footprint_objects.industry_characteristics = build_eora_industry_data(footprint_objects_params.eora_x_filename, footprint_objects_params.eora_ghg_filename, footprint_objects.un_to_iucn_codes); 
        elseif strcmp(system_type, 'hscpc')
            footprint_objects.industry_characteristics = build_hscpc_industry_data(footprint_objects_params.hscpc_x_filename, footprint_objects_params.hscpc_country_codes_filename, footprint_objects.un_to_iucn_codes);
        end
    
        footprint_objects.low_income_country_names = load_low_income_country_data(footprint_objects_params.low_income_countries_filename);
    
    elseif strcmp(input_type, 'footprint')
            
        footprint_objects.aquatic_data = partition_aquatic_data();
        marine_taxon_matches = footprint_objects.aquatic_data.marine_classified_taxon_matches;
        freshwater_taxon_matches = footprint_objects.aquatic_data.freshwater_classified_taxon_matches;
       
        footprint_objects.species_group_names = {'ANIMALIA', 'PLANTAE', 'marine', 'freshwater'}';
     
        footprint_objects.satellite_collapse_groups = cell(4, 1);
        footprint_objects.satellite_collapse_groups{1} = find(strcmp(satellite_species_characteristics.iucn_species_kingdom, footprint_objects.species_group_names{1}) ...
                                                                           & ~(marine_taxon_matches | freshwater_taxon_matches) );
        footprint_objects.satellite_collapse_groups{2} =  find(strcmp(satellite_species_characteristics.iucn_species_kingdom, footprint_objects.species_group_names{2}) ...
                                                                          & ~( marine_taxon_matches | freshwater_taxon_matches) ) ;
        footprint_objects.satellite_collapse_groups{3} = find(marine_taxon_matches);
        footprint_objects.satellite_collapse_groups{4} = find(freshwater_taxon_matches);
        
    end
        
end

function [un_to_iucn_codes] = load_un_to_iucn_codes(un_to_iucn_codes_filename)
    
    fid = fopen(un_to_iucn_codes_filename);
        un_to_iucn_data = textscan(fid,'%s %s %s', 'delimiter', ';');
    fclose(fid);
    un_to_iucn_codes = struct();
    un_to_iucn_codes.iucn_country_codes = un_to_iucn_data{3};
    un_to_iucn_codes.un_country_codes = un_to_iucn_data{2};
    un_to_iucn_codes.iucn_country_names = un_to_iucn_data{1};
    
end

function low_income_country_names = load_low_income_country_data(low_income_countries_filename)
    
    fid = fopen(low_income_countries_filename);
        low_income_country_names = textscan(fid,'%s', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);
    
end

function [industry_characteristics] = build_eora_industry_data(eora_x_filename, eora_ghg_filename, un_to_iucn_codes)       

    fid = fopen(eora_x_filename);
        x_data = textscan(fid,'%s %s %s %s %f', 'HeaderLines', 0, 'delimiter', ';');
    fclose(fid);
 
    [unique_countries, unique_inds] = unique(x_data{1}, 'stable');
    country_indexes = 1:length(unique_countries);
    country_index_list = zeros(size(x_data{1}));

    for country_ind = country_indexes
        country_index_list(strcmp(x_data{1}, unique_countries(country_ind))) = country_ind;
    end
    
    industry_characteristics = struct();
    industry_characteristics.country_names_list = x_data{1};
    industry_characteristics.country_codes_list = x_data{2};
    industry_characteristics.country_index_list = country_index_list;
    industry_characteristics.industry_commodity_discriminator_list = x_data{3};
    industry_characteristics.country_index_map = [unique_countries num2cell(country_indexes')];
    industry_characteristics.unique_countries = unique_countries;
    industry_characteristics.unique_country_codes = industry_characteristics.country_codes_list(unique_inds);
    
    industry_characteristics.commodity_classification_list = x_data{4};
    industry_characteristics.numerical_industry_data = x_data{5};
    industry_characteristics.industry_data_list = strcmp(industry_characteristics.industry_commodity_discriminator_list, 'Industries');
    industry_characteristics.commodity_data_list = strcmp(industry_characteristics.industry_commodity_discriminator_list, 'Commodities');
    
    [industry_characteristics.x_un, industry_characteristics.industry_codes_to_use] = build_iucn_industries_data(industry_characteristics.numerical_industry_data, ...
                                                                                                                 industry_characteristics.country_codes_list, ...
                                                                                                                 industry_characteristics.industry_data_list, ...
                                                                                                                 industry_characteristics.commodity_data_list, ...
                                                                                                                 un_to_iucn_codes);
                                    
    [industry_characteristics.x_un_names, ~] = build_iucn_industries_data(industry_characteristics.commodity_classification_list, ...
                                                                          industry_characteristics.country_codes_list, ...
                                                                          industry_characteristics.industry_data_list, ...
                                                                          industry_characteristics.commodity_data_list, ...
                                                                          un_to_iucn_codes); 
                                                                      
    industry_characteristics.ghg = build_eora_ghg(eora_ghg_filename, un_to_iucn_codes);
    
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


function [industry_characteristics] = build_hscpc_industry_data(hscpc_x_filename, hscpc_country_codes_filename, un_to_iucn_codes)
    
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


function ghg_un = build_eora_ghg(eora_ghg_filename, un_to_iucn_codes)
    
    fid = fopen(eora_ghg_filename);
        ghg_eora_data = textscan(fid,'%f %s %s %s %s %f', 'HeaderLines', 1, 'delimiter', ';');
    fclose(fid);
    
    industry_characteristics = struct();
    industry_characteristics.country_codes_list = ghg_eora_data{3};
    industry_data_list = strcmp(ghg_eora_data{4}, 'Industries');
    commodities_data_list = strcmp(ghg_eora_data{4}, 'Commodities');
    industry_characteristics.numerical_industry_data = ghg_eora_data{6};
    
    industry_names = ghg_eora_data{5};
    
%     export_data = regexpi(industry_names,'export');
%     non_export_data = cellfun('isempty', export_data);
    
    [ghg_un, ~] = build_iucn_industries_data(industry_characteristics.numerical_industry_data, industry_characteristics.country_codes_list, ...
                                                                industry_data_list, commodities_data_list, un_to_iucn_codes);
    
%     export_data = regexpi(industry_names,'export');
%     non_export_data = cellfun('isempty', export_data);
%     
%     
%     empty_industries = cellfun('isempty', ghg_industries);
%     empty_comms = cellfun('isempty', ghg_comms);
%     
%     inds_to_replace = empty_industries & ~empty_comms;
%     %data_to_use = industry_data & (ghg_CO2 > 0);
%     eora_ghg = ghg_industries;
%     eora_ghg(inds_to_replace) = ghg_comms(inds_to_replace);
    
end

   




function aquatic_identifiers = partition_aquatic_data()
    
    fid = fopen('~/GitHub/mrio_bio_satellite/iucn_input_data/improved_marine_identifiers.txt');
        aquatic_species_data = textscan(fid,'%s %f %s %s %f %s %f %s %f %s ', 'HeaderLines', 1, 'delimiter', '\t');
    fclose(fid);

    aquatic_data = aquatic_species_data{3};

    aquatic_identifiers = struct();
    [aquatic_identifiers.category_names, ~, aquatic_identifiers.categories] = unique( aquatic_data);
    aquatic_identifiers.frequencies = histc(aquatic_identifiers.categories, 1:numel(aquatic_identifiers.category_names));

    aquatic_identifiers.category_flags = cellfun(@(x) strcmp(aquatic_data, x), aquatic_identifiers.category_names, 'un', false);

    aquatic_identifiers.marine_ids = [5 6 8:10];
    aquatic_identifiers.freshwater_ids = [2:4 7];
    marine_classified_taxon_matches = aquatic_identifiers.category_flags([aquatic_identifiers.marine_ids]);
    freshwater_classified_taxon_matches = aquatic_identifiers.category_flags([aquatic_identifiers.freshwater_ids]);
    aquatic_identifiers.marine_classified_taxon_matches = sum([marine_classified_taxon_matches{:}], 2) > 0;
    aquatic_identifiers.freshwater_classified_taxon_matches = sum([freshwater_classified_taxon_matches{:}], 2) > 0;
    
end
