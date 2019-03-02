params_object = build_iucn_params();
processed_iucn_data_object = process_raw_iucn_data(params_object.satellite_params);
satellite_object = build_iucn_satellite(processed_iucn_data_object, params_object.satellite_params);
%% manfred script here 

load([params_object.process_MRIO_outputs_params.satellite_species_characteristics_filename]);
[~, ~, inds_to_use] = intersect(species_characteristics.species_taxons, iucn_data_object.iucn_threat_taxons, 'stable');
species_characteristics.species_kingdom = iucn_data_object.iucn_species_kingdom(inds_to_use);

trade_characteristics = analyse_global_consumption_routines(processed_iucn_data_object, analyse_mrio_params, species_characteristics);
