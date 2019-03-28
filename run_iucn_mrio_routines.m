iucn_params = build_iucn_params();

industry_inputs = build_input_objects('industry', iucn_params.footprint_objects_params, iucn_params.global_params.system_type);

processed_iucn_data = build_processed_iucn_data(iucn_params.satellite_params, ...
                                                industry_inputs.un_to_iucn_codes, ...
                                                industry_inputs.industry_characteristics, ...
                                                iucn_params.global_params);

satellite_object = build_iucn_satellite(processed_iucn_data, iucn_params.satellite_params, iucn_params.global_params.system_type, industry_inputs); clear processed_iucn_data

footprint_inputs = build_input_objects('footprint', iucn_params.footprint_objects_params, iucn_params.global_params.system_type, satellite_object.species_characteristics);

run_mrio_species_origin_destination(satellite_object.direct_satellite, footprint_inputs, industry_inputs, iucn_params.build_footprint_params, 'consumption', 2013);
                                                         
trade_characteristics = process_iucn_footprints(industry_inputs, ...
                                                footprint_inputs, ...
                                                iucn_params.analyse_footprint_params, ...
                                                iucn_params.build_footprint_params);  

                                            