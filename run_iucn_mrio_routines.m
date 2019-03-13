params_object = build_iucn_params();

processed_iucn_data = process_raw_iucn_data(params_object.satellite_params);
load([params_object.process_MRIO_outputs_params.iucn_data_object_filename]);
satellite_object = build_iucn_satellite(processed_iucn_data, params_object.satellite_params);
%% manfred script here 

load([params_object.process_MRIO_outputs_params.satellite_species_characteristics_filename]);
[~, ~, inds_to_use] = intersect(species_characteristics.species_taxons, iucn_data_object.iucn_threat_taxons, 'stable');
species_characteristics.species_kingdom = iucn_data_object.iucn_species_kingdom(inds_to_use);

display_satellite(satellite_object, params_object.satellite_params, species_characteristics, satellite_object.satellite_params.sector_lengths, satellite_object.satellite_params.sorted_country_names)

processed_finalsale_data = process_mrio_outputs_at_finalsale_level(processed_iucn_data, analyse_mrio_params, species_characteristics);
trade_characteristics = process_mrio_outputs_at_consumption_level(processed_iucn_data, analyse_mrio_params, species_characteristics, processed_finalsale_data);