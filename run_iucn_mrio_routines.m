iucn_params = build_iucn_params();

industry_inputs = build_input_objects('industry', iucn_params.footprint_objects_params, iucn_params.global_params.system_type);

satellite_inputs = build_processed_iucn_data(industry_inputs, ...
                                             iucn_params.satellite_params, ...
                                              iucn_params.global_params,  ...
                                             iucn_params.satellite_params);

satellite_object = build_iucn_satellite(satellite_inputs, iucn_params.satellite_params, iucn_params.global_params.system_type); 

footprint_inputs = build_input_objects('footprint', iucn_params.footprint_objects_params, iucn_params.global_params.system_type, satellite_inputs.species_characteristics);

run_mrio_species_origin_destination(satellite_object.direct_satellite, footprint_inputs, industry_inputs, iucn_params.build_footprint_params, 'consumption', 2013);
                                                         
footprints = process_iucn_footprints('global', 'consumption', industry_inputs, footprint_inputs, iucn_params.analyse_footprint_params, iucn_params.build_footprint_params);  

balanced_footprint = build_balanced_footprint(footprints, footprint_inputs, industry_inputs);

%writetable(footprints.ANIMALIA.aggregated_sector_scale.table, '~/Github/mrio_bio_satellite/footprint_outputs/Australia_consumption_aggregated_finalsale_sector_scale.txt', 'Delimiter', 'tab')

