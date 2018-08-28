%file_path = '~/Documents/MATLAB/BIO_SATELLITE/';  %%%%%% MODIFY THIS %%%%%%%%

function [IUCN_data_object, IUCN_params] = build_IUCN_data_object(file_path)
    
    tic
    addpath(file_path)
    IUCN_params = load_IUCN_params();
    IUCN_data_object = run_IUCN_data_routines(IUCN_params);
    disp(['IUCN data object built at ' num2str(round(100*toc/60)/100) ' mins, saving and writing tensors'])
    if ~exist(IUCN_params.Redlist_file_path, 'dir')
            mkdir(IUCN_params.Redlist_file_path); 
    end
    toc
    save([IUCN_params.Redlist_file_path, 'IUCN_data_object_2017.mat' ], 'IUCN_data_object');
    toc
    write_IUCN_causes(IUCN_data_object, IUCN_params.Redlist_file_path);
    %counters_check = false;
    % file_path = '~/Documents/MATLAB/BIO_SATELLITE/';  %%%%%% MODIFY THIS %%%%%%%%
    % addpath(file_path)
    % load('~/Documents/MATLAB/BIO_SATELLITE/IUCNRedList/IUCN_data_object_2016.mat')
    % load('~/Documents/MATLAB/BIO_SATELLITE/IUCNRedList/IUCN_tensor.mat')
    % satellite_params = build_satellite_params(IUCN_data_object);
    % [domestic_satellite, global_satellite, satellite_params] = build_IUCN_satellite(IUCN_data_object, IUCN_tensor, satellite_params);
end

