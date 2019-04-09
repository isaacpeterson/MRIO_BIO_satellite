function consumption_level_footprints = run_mrio_species_origin_destination(current_satellite, footprint_input_objects, industry_inputs, build_footprint_params, footprint_level, mrio_yr)
    
%     if ~build_footprint_params.build_footprint && exist(footprint_filename, 'file')
%         footprints = load(footprint_filename);
%         obj_names = fieldnames(footprints);
%         footprints = footprints.(obj_names{1});
%         return
%     end
    
    mrio_objects = build_mrio_objects(current_satellite, footprint_input_objects, build_footprint_params, mrio_yr, footprint_level);         

    disp(['Calculating footprints at ', footprint_level, ' level...'])
    
    if strcmp(footprint_level, 'finalsale')    
        
        mrio_objects.Y_to_use = sparse( sum(mrio_objects.Y, 2) ); 
        
        if strcmp(build_footprint_params.species_to_assess, 'all')
            species_to_assess = 1:size(mrio_objects.q, 1);
        else
            species_to_assess = build_footprint_params.species_to_assess;
        end
        
        for q_index = species_to_assess
            
            Ly = build_Ly(mrio_objects.L, mrio_objects.Y_to_use, build_footprint_params.use_sparse_representation);
            finalsale_footprint = build_current_footprint(mrio_objects.q(q_index, :), build_footprint_params, Ly, q_index);
            footprint_filename = [build_footprint_params.footprint_filename_prefix, footprint_level, '_', num2str(q_index), build_footprint_params.footprint_filename_suffix];
            save(footprint_filename, 'finalsale_footprint', '-v7.3')
            disp(['finalsale footprints for species ' num2str(q_index) ' done'])    
            
        end
        
    else
        
        if strcmp(build_footprint_params.countries_to_assess, 'all')
            countries_to_assess = 1:numel(industry_inputs.industry_characteristics.unique_countries);
        else
            [~, countries_to_assess] = intersect(industry_characteristics.unique_countries, build_footprint_params.countries_to_assess); 
        end
        
        inds_to_use = cellfun(@(x) ( (x - 1)*6 + 1): (x*6), num2cell(countries_to_assess), 'un', false);
        mrio_objects.Y_to_use = cellfun(@(x) sum( mrio_objects.Y(:, x), 2) + sum(mrio_objects.Y(:, (mrio_objects.y_length + x)), 2), inds_to_use, 'un', false);
        
        if ~exist(build_footprint_params.footprint_filename_prefix, 'dir')
            mkdir(build_footprint_params.footprint_filename_prefix)
        end
        
        for country_ind = 1:numel(countries_to_assess)
            
            Ly = build_Ly(mrio_objects.L, mrio_objects.Y_to_use{country_ind}, build_footprint_params.use_sparse_representation);
            consumption_level_footprints = struct();
            [consumption_level_footprints.production, consumption_level_footprints.finalsale, consumption_level_footprints.species, consumption_level_footprints.vals] = ...
                cellfun(@(x) build_current_footprint(mrio_objects.q(x, :), build_footprint_params, Ly, x), num2cell((1:length(footprint_input_objects.satellite_collapse_groups))'), 'un', false);
            
            save([build_footprint_params.footprint_filename_prefix footprint_level '_' industry_inputs.industry_characteristics.unique_country_codes{countries_to_assess(country_ind)}  ...
                  build_footprint_params.footprint_filename_suffix], 'consumption_level_footprints', '-v7.3')
            
            disp([industry_inputs.industry_characteristics.unique_country_codes{countries_to_assess(country_ind)}, ' consumption level footprints done'])  
            
        end
        
        
        
    end
    
end



function mrio_objects = build_mrio_objects(current_satellite, footprint_input_objects, build_footprint_params, mrio_yr, footprint_level)
        
    if (build_footprint_params.build_mrio_outputs == false) && exist([build_footprint_params.mrio_objects_file_prefix num2str(mrio_yr) '.mat'], 'file')
        disp(['loading processed MRIO data for year ' num2str(mrio_yr) '...']);
        load([build_footprint_params.mrio_objects_file_prefix num2str(mrio_yr) '.mat']); 
    else
        
        disp(['Processing MRIO data for year ' num2str(mrio_yr)]);
                
        mrio_objects = load_mrio_objects(build_footprint_params, mrio_yr);
                
        % 2.2.1 Cut off Statistical Discrepancies
       	mrio_objects = trim_mrio_objects(mrio_objects);
        mrio_objects.y_length = size(mrio_objects.Y, 2);        
        v_split = mirror_mrio( mrio_objects.V );
        y_split = mirror_mrio( mrio_objects.Y );
                
        mrio_objects.V = vertcat(v_split.pos, y_split.neg); 
        mrio_objects.Y = horzcat(y_split.pos, v_split.neg);
        
        mrio_objects.x = sum(mrio_objects.T, 2) + sum(mrio_objects.Y,2);
        
        disp('building A')
        mrio_objects.A = mrio_objects.T * diag( 1 ./ (mrio_objects.x + 1e-10) ); 
        
        disp('building Leontief inverse')

        mrio_objects.L = inv(eye(size(mrio_objects.A)) - mrio_objects.A);  
      
        save([build_footprint_params.mrio_objects_file_prefix num2str(mrio_yr) '.mat'], 'mrio_objects', '-v7.3');  
    	
    end
    
    if ~strcmp(footprint_level, 'finalsale')
    	current_satellite = run_collapse_satellite(current_satellite, footprint_input_objects.satellite_collapse_groups);
    end
        
    mrio_objects.q = current_satellite * diag(sparse(1 ./ (mrio_objects.x + 1e-10) ));
        
	if build_footprint_params.use_sparse_representation 
        	mrio_objects.q = sparse(mrio_objects.q); 
    end
        
         
end


function collapsed_satellite = run_collapse_satellite(current_satellite, satellite_collapse_characteristics)
    
    collapsed_satellite = cellfun(@(x) sum(current_satellite(x, :), 1), satellite_collapse_characteristics, 'un', false);
    collapsed_satellite = vertcat( collapsed_satellite{:});
    
end

function current_object = read_obj_from_structure(obj_filename)
    current_object = load(obj_filename);
    obj_names = fieldnames(current_object);
    current_object = current_object.(obj_names{1});
end

function Ly = build_Ly(L, Y, use_sparse_representation)
    
    if use_sparse_representation
    	Ly = sparse(L * sparse(diag(Y))); 
    else Ly = L * sparse(diag(Y)); 
    end
    
end


function [footprints_production, footprints_finalsale, footprints_species, footprints_vals] = build_current_footprint(current_qrow, build_footprint_params, Ly, q_index)

    current_footprint = sparse( diag(current_qrow) ) * Ly;
    
    if (build_footprint_params.thresh > 0)
        current_footprint = current_footprint.*(current_footprint > build_footprint_params.thresh);  
    end
    
    [footprints_production, footprints_finalsale, footprints_vals] = find(current_footprint); 
    
    footprints_species = repmat(q_index, [size(footprints_vals, 1) 1]);
        
end


function mrio_objects = load_mrio_objects(build_footprint_params, mrio_yr)

    mrio_objects = struct();

    mrio_objects.Y = binread([build_footprint_params.mrio_data_file_prefix '_Y-Results_' num2str(mrio_yr), build_footprint_params.mrio_data_file_suffix]);
    mrio_objects.T = binread([build_footprint_params.mrio_data_file_prefix '_T-Results_' num2str(mrio_yr), build_footprint_params.mrio_data_file_suffix]);
    mrio_objects.V = double(full(binread([build_footprint_params.mrio_data_file_prefix '_V-Results_' num2str(mrio_yr), build_footprint_params.mrio_data_file_suffix]))); 
    
end

function mrio_objects = trim_mrio_objects(mrio_objects)
    
    trim_vec = 1:(size(mrio_objects.T, 1) - 1);
    mrio_objects.V = mrio_objects.V(:, trim_vec);
    mrio_objects.Y = mrio_objects.Y(trim_vec, :);
    mrio_objects.T = mrio_objects.T( trim_vec, trim_vec);
        
end            

function v_mirror = mirror_mrio(V)
    
    v_mirror = struct();
    v_mirror.pos = V .* (V > 0);
    v_mirror.neg = -(V .* (V < 0))';
    
end


function footprint = run_exclusions(footprint)
    
    disp('Processing exclusions.');
        % Origin country exclusions  
      
    exclusions.production_country_exclusions = ismember(footprint.production, find(ismember(build_footprint_params.eora_full_labels(:, 1), build_footprint_params.production_country_exclusions)));
    exclusions.production_industry_exclusions = ismember(footprint.production, find(ismember(build_footprint_params.eora_full_labels(:, 2), build_footprint_params.production_industry_exclusions)));
    exclusions.finalsale_country_exclusions = ismember(footprint.finalsale, find(ismember(build_footprint_params.eora_full_labels(:, 1), build_footprint_params.finalsale_country_exclusions)));
    exclusions.finalsale_industry_exclusions = ismember(footprint.finalsale, find(ismember(build_footprint_params.eora_full_labels(:, 2), build_footprint_params.finalsale_industry_exclusions)));

    footprint = structfun(x(~exclusions), footprint, exclusions);
    
end


% function current_mrio_object = load_current_mrio_object(current_filename, backup_filename)
% 
%     if exist(current_filename,'file');
%     	current_mrio_object = binread(current_filename);
%     else
%         current_mrio_object = binread(backup_filename);
%     end;
%     
% end
% 
% function mrio_labels = build_mrio_labels()
% 
%     disp(['Loading and processing labels and acronyms.']);
% 
%     % 2.1 Load tourism final demand
%     % 2.1.1 Load country names
%     if exist([concdir 'countryNames.mat'],'file');
%         load([concdir 'countryNames.mat']);
%     else
%         [~,countryNames] = xlsread([build_footprint_params.mrio_data_dir 'Tourism-Z_' num2str(2013) '.xlsx'],'m','A2:A190');
%         save([concdir 'countryNames.mat'],'countryNames');
%     end;
%     
%     NCOUN = size(countryNames,1);
%     % 2.1.2 Load country acronyms
%     if exist([concdir 'countryAcros.mat'],'file');
%         load([concdir 'countryAcros.mat']);
%     else
%         [~,countryAcros] = xlsread([build_footprint_params.mrio_data_dir 'SumOfWTOdata.xlsx'],'Pop','B2:B190');
%         save([concdir 'countryAcros.mat'],'countryAcros');
%     end;
%     
%     countryAcros{27,1} = 'VGB'; % correct British Virgin Islands
%     % 2.1.4 Load labels
%     
%     load([concdir 'EoraFullLabels.mat']); load([concdir 'Eora26Labels.mat']);
%     
%     if exist([concdir 'EoraSectorNumbers.mat'],'file');
%         load([concdir 'EoraSectorNumbers.mat']);
%     else
%         [~,~,EoraSectorNumbers] = xlsread([treedir 'EoraSectorNum.xlsx'],'EoraSectorNum');
%         save([concdir 'EoraSectorNumbers.mat'],'EoraSectorNumbers');
%     end; 
%     
%     if exist([concdir 'EoraFullAcronymsEntities.mat'],'file');
%         load([concdir 'EoraFullAcronymsEntities.mat']);
%     else
%         [~,EoraFullAcronymsEntities] = xlsread([treedir 'index_t.xlsx'],'index_t','C2:D14839');
%         save([concdir 'EoraFullAcronymsEntities.mat'],'EoraFullAcronymsEntities');
%     end; 
%     
%     load([workdir 'satellite_species_params.mat']);
%     % SpecName = cellfun(@(x,y,z) [ x ' (' y '; ' z ')' ] , sat_params.species_names, sat_params.species_classification, num2cell(sat_params.species_taxons), 'UniformOutput', false);
%     SpecName = cellfun(@(x,y) [ x ' (' y ')' ] , sat_params.species_names, sat_params.species_classification, 'UniformOutput', false);
% 
% end
% 
