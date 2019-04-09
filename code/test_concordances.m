load('~/GitHub/MRIO_BIO_satellite/concordances/eora_threat_concordances/eora_threat_concordance_group.mat')
load('~/GitHub/MRIO_BIO_satellite/processed_data/eora/eora_satellite_inputs.mat')
load('/Users/isaacpet/GitHub/MRIO_BIO_satellite/iucn_input_data/EoraFullLabels.mat')

fid = fopen('~/GitHub/MRIO_BIO_satellite/iucn_input_data/threatcauseclassification.txt');
	old_threat_cause_class_data = textscan(fid,'%s %s %f', 'HeaderLines', 0, 'delimiter', ';');
fclose(fid);

iucn_threat_type = strrep(old_threat_cause_class_data{1}, "'", "");

fid = fopen('~/GitHub/MRIO_BIO_satellite/iucn_input_data/x_data_187.txt');
    industry_data = textscan(fid,'%s %s %s %s %f', 'HeaderLines', 0, 'delimiter', ';');
fclose(fid);
    
write_to_table = false;
overwrite_eora26 = true;

lengths = cellfun(@(x) size(x, 2), threat_concordance);
countries_to_replace = find(lengths == 25);

build_barplot = false;
overwrite_barplot = true;

if build_barplot == true
    %data_to_use = ismember(processed_data.iucn_data.iucn_species_data.country_names_list, satellite_inputs.country_characteristics.country_names(countries_to_replace));
    data_to_use = 1:numel(processed_data.iucn_data.iucn_species_data.country_names_list);


    [unique_names, ~, J] = unique(processed_data.iucn_data.iucn_species_data.iucn_threats_list(data_to_use));
    unique_counts_balanced = histc(J, 1:numel(unique_names));


    [~, b, c] = intersect(iucn_threat_type, unique_names);

    names_to_replace = unique_names;
    names_to_replace(c) = old_threat_cause_class_data{2}(b);


    iucn_threat_type_to_use = setdiff(1:length(iucn_threat_type), 9:84); 

    [unique_names, ~, J] = unique(processed_data.iucn_data.iucn_threat_data.iucn_threats);
    unique_counts_unbalanced = histc(J, 1:numel(unique_names));

    FigH = figure('Position', get(0, 'Screensize'));

    bar(unique_counts_balanced/sum(unique_counts_balanced), 'FaceColor', 'red')

    if overwrite_barplot == true
        hold on
        bar(unique_counts_unbalanced/sum(unique_counts_unbalanced), 0.5, 'FaceColor', 'blue')
        hold off
    end
    set(gca, 'xtick', 1:numel(unique_names))
    set(gca, 'XTickLabel', names_to_replace)
    set(gca, 'XTickLabelRotation', 90)
    F = getframe(FigH);
    imwrite(F.cdata, '~/GitHub/MRIO_BIO_satellite/global_threat_histogram_balanced_unbalanced.png', 'png')
    close(FigH)

end

if overwrite_eora26 == true
     
    countries_to_write = countries_to_replace';
    industry_names = matlab.lang.makeValidName(satellite_inputs.country_characteristics.x_names{1});

    fid = fopen('~/GitHub/MRIO_BIO_satellite/concordances/eora_threat_concordances/human_readable/eora_26_template_corrected.txt');
        eora_26_template = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1, 'delimiter', 'tab');
    fclose(fid);
    current_concordance = [eora_26_template{:}]; 
else 
    countries_to_write = setdiff(1:numel(threat_concordance), countries_to_replace);
end


for country_ind = countries_to_write

    current_concordance_filename = ['~/GitHub/MRIO_BIO_satellite/concordances/eora_threat_concordances/human_readable/concordances_to_correct/', ...
                                    satellite_inputs.country_characteristics.country_names{country_ind}, '.txt'];
                
    if write_to_table == true
        
        if ~overwrite_eora26
            current_concordance = threat_concordance{country_ind};
        end
        
        industry_names = matlab.lang.makeValidName(satellite_inputs.country_characteristics.x_names{country_ind});
        
        [unique_names, ~, J] = unique(industry_names);
        unique_counts_balanced = find(histc(J, 1:numel(unique_names)) > 1);
        
        for replace_ind = 1:numel(unique_counts_balanced)
            duplicate_inds = find(strcmp(industry_names, unique_names(unique_counts_balanced(replace_ind))));
            industry_names{duplicate_inds(end)} = strcat([industry_names{duplicate_inds(end)} 'b']);
        end
        
        T = cell2table([old_threat_cause_class_data{2} num2cell(current_concordance)], 'VariableNames', [{'IUCN_threat_type'} industry_names(1:(end - 1))']);
        
        writetable(T, current_concordance_filename, 'delimiter', 'tab');
         
    else
        
        current_concordance = dlmread(current_concordance_filename, '\t', 1, 1);              

        csvwrite(['~/GitHub/MRIO_BIO_satellite/concordances/eora_threat_concordances/human_readable/corrected_concordances/', ...
                    satellite_inputs.country_characteristics.country_codes{country_ind}, '.csv'], current_concordance);                
    end
    
end