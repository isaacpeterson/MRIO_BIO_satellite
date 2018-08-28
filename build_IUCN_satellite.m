IUCN_data_params = load_IUCN_params();
IUCN_data_object = run_IUCN_data_routines(IUCN_data_params);

satellite_params = build_satellite_params(IUCN_data_object);

[domestic_satellite, global_satellite, satellite_params] = run_IUCN_satellite_routines_HSCPC(IUCN_data_object, IUCN_data_params, satellite_params);

%save('~/Documents/MATLAB/BIO_SATELLITE/RedList_2016/total_domestic_satellite.mat', 'domestic_satellite', '-v7.3')
%save('~/Documents/MATLAB/BIO_SATELLITE/RedList_2016/total_global_satellite.mat', 'global_satellite', '-v7.3')
%save('~/Documents/MATLAB/BIO_SATELLITE/RedList_2016/satellite_params.mat', 'global_satellite', '-v7.3')

% display_satellite(domestic_satellite, global_satellite, satellite_params)
% 
% au_inds = find(strcmp(IUCN_data_object.IUCN_country_codes_list, 'AU'));
% [a, b] = unique(IUCN_data_object.IUCN_taxons_list(au_inds));
% 
% IUCN_data_object.IUCN_status_inds(au_inds(b))



