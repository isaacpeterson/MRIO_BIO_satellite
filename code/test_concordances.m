load('~/GitHub/MRIO_BIO_satellite/concordances/eora_threat_concordances/eora_threat_concordance_group.mat')
load('~/GitHub/MRIO_BIO_satellite/processed_data/eora/eora_satellite_inputs.mat')

fid = fopen('~/GitHub/MRIO_BIO_satellite/iucn_input_data/threatcauseclassification.txt');
	old_threat_cause_class_data = textscan(fid,'%s %s %f', 'HeaderLines', 0, 'delimiter', ';');
fclose(fid);

bad_chars = {'&', ' ', '&', '-', ','};

for country_ind = 1:numel(threat_concordance)
    industry_names = matlab.lang.makeValidName(satellite_inputs.iucn_industry_characteristics.x_names{country_ind});

%     for bad_ind = 1:length(bad_chars)
%         industry_names = strrep(industry_names, bad_chars{bad_ind}, "");
%     end
    
    cell2table([old_threat_cause_class_data{1} num2cell(threat_concordance{country_ind})], 'VariableNames', industry_names)
end