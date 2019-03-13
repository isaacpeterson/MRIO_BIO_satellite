function run_mrio_species_origin_destination(mrio_footprint_params)

    %mrio_labels = build_mrio_labels();
    
    if strcmp(mrio_footprint_params.assessment_scale, 'finalsale')
        country_vec = 1;
    else
        country_vec = mrio_footprint_params.country_vec;
    end
    
    for mrio_yr = mrio_footprint_params.startyear:mrio_footprint_params.endyear 
        
        mrio_objects = build_mrio_objects(mrio_footprint_params, mrio_yr);
        
        disp('Calculating disaggregated footprints...')
    
        for consumption_country_ind = country_vec
            build_footprints(mrio_objects, mrio_footprint_params, consumption_country_ind);
        end
    
    end
    
end



function mrio_objects = build_mrio_objects(mrio_footprint_params, mrio_yr)
        
    if (mrio_footprint_params.build_mrio_outputs == false) && exist([mrio_footprint_params.rawdatadir 'mrio_objects_' num2str(mrio_yr) '.mat'], 'file')
        disp(['loading financial data for year ' num2str(mrio_yr) '...']);
        load([mrio_footprint_params.rawdatadir 'mrio_objects_' num2str(mrio_yr) '.mat']); 
    else
        
        disp(['Processing raw MRIO data for year ' num2str(mrio_yr) '...']);
                
        mrio_objects = load_mrio_objects(mrio_footprint_params, mrio_yr);
                
        % 2.2.1 Cut off Statistical Discrepancies
       	mrio_objects = trim_mrio_objects(mrio_objects);
                
        v_split = mirror_mrio( mrio_objects.V );
        y_split = mirror_mrio( mrio_objects.Y );
                
        mrio_objects.V = vertcat(v_split.pos, y_split.neg); 
        mrio_objects.Y = horzcat(y_split.pos, v_split.neg);
        mrio_objects.x = sum(mrio_objects.T, 2) + sum(mrio_objects.Y,2);
        
        disp('calculating A')
        mrio_objects.A = mrio_objects.T * diag( 1 ./ (mrio_objects.x + 1e-10) ); 
        
        disp('calculating Leontief inverse')
        % 2.2.3 Calculate Leontief inverse  
        mrio_objects.L = inv(eye(size(mrio_objects.A)) - mrio_objects.A);  
        
        current_satellite = read_obj(mrio_footprint_params.satellite_filename);
        if mrio_footprint_params.use_sparse_representation 
        	mrio_objects.q = sparse(current_satellite * diag( sparse( 1./ (mrio_objects.x + 1e-10) ) ) ); 
        else mrio_objects.q = current_satellite * diag(sparse(1 ./ (mrio_objects.x + 1e-10) ) ); 
        end;
                
        save([mrio_footprint_params.rawdatadir 'mrio_objects_' num2str(mrio_yr) '.mat'], 'mrio_objects', '-v7.3');  
    	
    end
    
	if strcmp(mrio_footprint_params.assessment_scale, 'finalsale')
    	mrio_objects.Y_to_use = sparse( sum(mrio_objects.Y, 2) );
    else
        mrio_objects.Y_to_use = cellfun(@(x) sum( mrio_objects.Y(:, (x -1)*6 + 1:x*6), 2) + sum(mrio_objects.Y(:, 1140 + [ (x-1)*6+1:x*6] ), 2), num2cell(mrio_footprint_params.q_vec));
    end
              
end

function current_object = read_obj(obj_filename)
    current_object = load(obj_filename);
    obj_names = fieldnames(current_object);
    current_object = current_object.(obj_names{1});
end

function Ly = build_Ly(L, Y, mrio_footprint_params)

    if mrio_footprint_params.use_sparse_representation; 
    	Ly = sparse(L * sparse(diag(Y))); 
    else Ly = L * sparse(diag(Y)); 
    end;
    
end


function footprint = build_footprints(mrio_objects, mrio_footprint_params, consumption_country_ind)

	if strcmp(mrio_footprint_params.assessment_scale, 'finalsale')
        Ly = build_Ly(mrio_objects.L, mrio_objects.Y_to_use, mrio_footprint_params);
    else
    	Ly = build_Ly(mrio_objects.L, mrio_objects.Y_to_use{consumption_country_ind}, mrio_footprint_params);
    end
    
	footprint = struct();
	[footprint.production, footprint.finalsale, footprint.vals, footprint.species] = cellfun(@(x) build_current_footprint(mrio_objects.q(mrio_footprint_params.q_vec(x), :), mrio_footprint_params, Ly, x), ...
                                                                                             num2cell(1:numel(mrio_footprint_params.q_vec)), 'un', false);
	
    footprint = structfun(@(x) vertcat(x{:}), footprint, 'un', false);
        
	if strcmp(mrio_footprint_params.assessment_scale, 'finalsale')
        footprint_filename = [mrio_footprint_params.workdir mrio_footprint_params.satellite_type '_footprint_finalsale.mat'];
    else
        footprint.consumption_country = repmat(consumption_country_ind, [numel(footprint.vals) 1]);
        footprint_filename = [mrio_footprint_params.workdir mrio_footprint_params.country_names{consumption_country_ind} '_' mrio_footprint_params.satellite_type '_footprint_consumption_level.mat'];
    end
    
    save(footprint_filename, 'footprint');
    
end

function [footprint_production, footprint_finalsale, footprint_vals, footprint_species] = build_current_footprint(current_qrow, mrio_footprint_params, Ly, q_index)

    current_footprint = sparse( diag(current_qrow) ) * Ly;
    current_footprint = current_footprint .* (current_footprint > mrio_footprint_params.thresh);               
    [footprint_production, footprint_finalsale, footprint_vals] = find(current_footprint);   
    footprint_species = repmat(q_index, [size(footprint_vals, 1) 1]);
    
end


function mrio_objects = load_mrio_objects(mrio_footprint_params, mrio_yr)

    mrio_objects = struct();

    mrio_objects.Y = binread([mrio_footprint_params.mrio_outputs_file_prefix '_Y-Results_' num2str(mrio_yr), mrio_footprint_params.mrio_outputs_file_suffix]);
    mrio_objects.T = binread([mrio_footprint_params.mrio_outputs_file_prefix '_T-Results_' num2str(mrio_yr), mrio_footprint_params.mrio_outputs_file_suffix]);
    mrio_objects.V = double(full(binread([mrio_footprint_params.mrio_outputs_file_prefix '_V-Results_' num2str(mrio_yr), mrio_footprint_params.mrio_outputs_file_suffix]))); 
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
    v_neg = V .* (V < 0);
    v_mirror.neg = -(v_neg)';
    
end


function footprint = run_exclusions(footprint)
    
    disp('Processing exclusions.');
        % Origin country exclusions  
      
    exclusions.production_country_exclusions = ismember(footprint.production, find(ismember(mrio_footprint_params.eora_full_labels(:, 1), mrio_footprint_params.production_country_exclusions)));
    exclusions.production_industry_exclusions = ismember(footprint.production, find(ismember(mrio_footprint_params.eora_full_labels(:, 2), mrio_footprint_params.production_industry_exclusions)));
    exclusions.finalsale_country_exclusions = ismember(footprint.finalsale, find(ismember(mrio_footprint_params.eora_full_labels(:, 1), mrio_footprint_params.finalsale_country_exclusions)));
    exclusions.finalsale_industry_exclusions = ismember(footprint.finalsale, find(ismember(mrio_footprint_params.eora_full_labels(:, 2), mrio_footprint_params.finalsale_industry_exclusions)));

    footprint = structfun(x(~exclusions), footprint, exclusions);
    
end


function current_mrio_object = load_current_mrio_object(current_filename, backup_filename)

    if exist(current_filename,'file');
    	current_mrio_object = binread(current_filename);
    else
        current_mrio_object = binread(backup_filename);
    end;
    
end

function mrio_labels = build_mrio_labels()

    disp(['Loading and processing labels and acronyms.']);

    % 2.1 Load tourism final demand
    % 2.1.1 Load country names
    if exist([concdir 'countryNames.mat'],'file');
        load([concdir 'countryNames.mat']);
    else
        [~,countryNames] = xlsread([mrio_footprint_params.rawdatadir 'Tourism-Z_' num2str(2013) '.xlsx'],'m','A2:A190');
        save([concdir 'countryNames.mat'],'countryNames');
    end;
    
    NCOUN = size(countryNames,1);
    % 2.1.2 Load country acronyms
    if exist([concdir 'countryAcros.mat'],'file');
        load([concdir 'countryAcros.mat']);
    else
        [~,countryAcros] = xlsread([mrio_footprint_params.rawdatadir 'SumOfWTOdata.xlsx'],'Pop','B2:B190');
        save([concdir 'countryAcros.mat'],'countryAcros');
    end;
    countryAcros{27,1} = 'VGB'; % correct British Virgin Islands
    % 2.1.4 Load labels
    
    load([concdir 'EoraFullLabels.mat']); load([concdir 'Eora26Labels.mat']);
    
    if exist([concdir 'EoraSectorNumbers.mat'],'file');
        load([concdir 'EoraSectorNumbers.mat']);
    else
        [~,~,EoraSectorNumbers] = xlsread([treedir 'EoraSectorNum.xlsx'],'EoraSectorNum');
        save([concdir 'EoraSectorNumbers.mat'],'EoraSectorNumbers');
    end; 
    
    if exist([concdir 'EoraFullAcronymsEntities.mat'],'file');
        load([concdir 'EoraFullAcronymsEntities.mat']);
    else
        [~,EoraFullAcronymsEntities] = xlsread([treedir 'index_t.xlsx'],'index_t','C2:D14839');
        save([concdir 'EoraFullAcronymsEntities.mat'],'EoraFullAcronymsEntities');
    end; 
    
    load([workdir 'satellite_species_params.mat']);
    % SpecName = cellfun(@(x,y,z) [ x ' (' y '; ' z ')' ] , sat_params.species_names, sat_params.species_classification, num2cell(sat_params.species_taxons), 'UniformOutput', false);
    SpecName = cellfun(@(x,y) [ x ' (' y ')' ] , sat_params.species_names, sat_params.species_classification, 'UniformOutput', false);

end

