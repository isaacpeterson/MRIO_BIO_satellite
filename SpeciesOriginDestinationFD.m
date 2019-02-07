%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% This code post- processes the updated species %%
%% threat satellite and issues results for the 
%% Isaac and Matthew
%%
%% Written by Manfred Lenzen
%% 7. Maerz 2018
%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1 Initialisation

clear; clc; close all;
if strcmp(computer,'MACI64');
elseif strcmp(computer,'GLNXA64');
    workdir = '/import/emily1/isa/Projects/2018_IUCN_update/';
    rawdatadir = ['/import/emily1/isa/Projects/2017_Tourism/FutuTourism/RawData/']; treedir = [rawdatadir 'TreeStructure/'];
    eoradir = '/import/laika1/isa/Eora/Phase199/Loop082/Results/';
    concdir = ['/import/emily1/isa/Projects/2017_Tourism/FutuTourism/Concordances/'];
end;
STARTYEAR = 2013; ENDYEAR = 2013; SP = 0;
scrsz = get(groot,'ScreenSize');
THRESH(1,1) = 0.001; THRESH(2,1) = 0.01;  
satstr{1,1} = 'global'; satstr{2,1} = 'domestic';
TOPNUMBER = 1000;

for sat = 1:2; % Loop over global and domestic satellites

    %% 2 Load data
    disp(['Loading and processing labels and acronyms.']);

    % 2.1 Load tourism final demand
    % 2.1.1 Load country names
    if exist([concdir 'countryNames.mat'],'file');
        load([concdir 'countryNames.mat']);
    else
        [~,countryNames] = xlsread([rawdatadir 'Tourism-Z_' num2str(2013) '.xlsx'],'m','A2:A190');
        save([concdir 'countryNames.mat'],'countryNames');
    end;
    NCOUN = size(countryNames,1);
    % 2.1.2 Load country acronyms
    if exist([concdir 'countryAcros.mat'],'file');
        load([concdir 'countryAcros.mat']);
    else
        [~,countryAcros] = xlsread([rawdatadir 'SumOfWTOdata.xlsx'],'Pop','B2:B190');
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
    load('species_kingdoms.mat');
    % SpecName = cellfun(@(x,y,z) [ x ' (' y '; ' z ')' ] , sat_params.species_names, sat_params.species_classification, num2cell(sat_params.species_taxons), 'UniformOutput', false);
    % SpecName = cellfun(@(x,y) [ x ' (' y ')' ] , sat_params.species_names, sat_params.species_classification, 'UniformOutput', false);
    SpecName = unique(species_kingdoms)'; % {'Erste';'Zweite';'Dritte';'Vierte'};
    
    if ~exist([workdir 'SpThrVals_domestic_1st4_Qrows.mat'],'file');
        
        % 2.2 Load IO matrices
        disp(['Loading and processing raw monetary data.']);
        for yy = STARTYEAR : ENDYEAR; % STARTYEAR : ENDYEAR;
            if exist([rawdatadir 'L_' num2str(yy) '.mat'],'file');
                load([rawdatadir 'L_' num2str(yy) '.mat']); 
                load([rawdatadir 'x_' num2str(yy) '.mat']); 
                load([rawdatadir 'Y_' num2str(yy) '.mat']);
            else
                disp(['Processing input-output data for year ' num2str(yy) '.']);
                if exist([eoradir '21160000_annetest_AllCountries_199_T-Results_' num2str(yy) '_082_Markup001.spbin'],'file');
                    ID = binread([eoradir '21160000_annetest_AllCountries_199_T-Results_' num2str(yy) '_082_Markup001.spbin']);
                else
                    ID = binread([eoradir '21140000_annetest_AllCountries_199_T-Results_' num2str(yy) '_082_Markup001.spbin']);
                end;
                if exist([eoradir '21160000_annetest_AllCountries_199_Y-Results_' num2str(yy) '_082_Markup001.spbin'],'file');
                    FD = binread([eoradir '21160000_annetest_AllCountries_199_Y-Results_' num2str(yy) '_082_Markup001.spbin']);
                else
                    FD = binread([eoradir '21140000_annetest_AllCountries_199_Y-Results_' num2str(yy) '_082_Markup001.spbin']);
                end;

                if exist([eoradir '21160000_annetest_AllCountries_199_V-Results_' num2str(yy) '_082_Markup001.spbin'],'file');
                    VA = binread([eoradir '21160000_annetest_AllCountries_199_V-Results_' num2str(yy) '_082_Markup001.spbin']);
                else
                    VA = binread([eoradir '21140000_annetest_AllCountries_199_V-Results_' num2str(yy) '_082_Markup001.spbin']);
                end;
                % 2.2.1 Cut off Statistical Discrepancies
                VA = VA(:,1:size(ID,1)-1); FD = FD(1:size(ID,1)-1,:); ID = ID(1:size(ID,1)-1,1:size(ID,1)-1); 
                % 2.2.2 Mirror
                VAneg = VA; VAneg(VAneg>0) = 0; FDneg = FD; FDneg(FDneg>0) = 0;
                VA(VA<0) = 0; FD(FD<0) = 0;
                VA = [VA ; -FDneg']; FD = [FD , -VAneg'];
                TO = sum(ID,2) + sum(FD,2); 
                % 2.2.3 Calculate Leontief inverse
                A = ID * diag(1./(TO+1e-10)); 
                L = inv(eye(size(A)) - A);       
                save([rawdatadir 'Y_' num2str(yy) '.mat'],'FD');  clear FD;
                save([rawdatadir 'T_' num2str(yy) '.mat'],'ID');  clear ID;
                save([rawdatadir 'A_' num2str(yy) '.mat'],'A');  clear A;
                save([rawdatadir 'x_' num2str(yy) '.mat'],'TO'); 
                save([rawdatadir 'L_' num2str(yy) '.mat'],'L');   
                clear FD VA FDneg VAneg;
            end;
        end;

        % 2.3 Compress final demand
        for i=1:190;
            tmpFD(:,i) = sum(FD(:,(i-1)*6+1:i*6),2) + sum(FD(:,1140+[(i-1)*6+1:i*6]),2);
        end;
        FD = full(tmpFD); clear tmpFD;   
        %FD = sparse(sum(FD,2));

        % 2.4 Load species threat data
        disp(['Loading and processing species threat data.']);
        % SpThr = rand([41000,size(L,1)]);
        if sat == 1;
            load([workdir 'total_global_satellite.mat']); SpThr = global_satellite; clear global_satellite;
        elseif sat == 2;
            load([workdir 'total_domestic_satellite.mat']); SpThr = domestic_satellite; clear domestic_satellite;
        end;
        if SP; q = sparse(SpThr * diag(sparse(1./(TO+1e-10)))); else q = SpThr * diag(sparse(1./(TO+1e-10))); end; clear SpThr;

        % 2.5 Calculate disabggregated footprints q # L # y
        % if exist(,'file');
        % else
            disp(['Calculating disaggregated footprints.']);
            SpThrSubs = []; SpThrVals = []; SpThrCountries = [];
            for sp = 1:size(SpecName,2);
                SpIndices = strcmp(species_kingdoms, SpecName{1,sp});
                qtmp(sp,:) = sum(q(SpIndices,:));
            end;
            q = qtmp; clear sp qtmp;
 
            for spthr = 1:size(q,1);
                for c=1:189;
                    disp(['Q row ' num2str(spthr) ', ' SpecName(spthr) '; country ' num2str(c)]);
                    FDc = FD(:,c);
                    if SP; Ly = sparse(L * sparse(diag(FDc))); else Ly = L * sparse(diag(FDc)); end;
                    tic;
                    tmp = q(spthr,:);
                    tmp = sparse(diag(tmp)) * Ly;
                    tmp(tmp<THRESH(sat,1)) = 0;
                    % Ftensor = sptensor(tmp); tmps = Ftensor.subs; tmpv = Ftensor.vals;
                    [tmpr tmpc,tmpv] = find(tmp);
                    SpThrSubs = [SpThrSubs ; uint16([repmat(spthr,[size(tmpv,1) 1]), tmpr, tmpc])]; 
                    SpThrVals = [SpThrVals; single(tmpv)];
                    SpThrCountries = [SpThrCountries; uint16(repmat(c,[size(tmpr,1) 1]))];
                end;
                %if floor(spthr/1000) == spthr/1000;
                    save([workdir 'SpThrSubs_' satstr{sat,1} '_1st' num2str(spthr) '_Qrows.mat'],'SpThrSubs');
                    save([workdir 'SpThrVals_' satstr{sat,1} '_1st' num2str(spthr) '_Qrows.mat'],'SpThrVals');
                    save([workdir 'SpThrCnts_' satstr{sat,1} '_1st' num2str(spthr) '_Qrows.mat'],'SpThrCountries');
                %end
                disp(['Threat satellite row ' num2str(spthr) ' processed in ' num2str(toc) ' s; found ' num2str(size(tmpv,1)) ' values.']); 
            end;

    else
        
        load([workdir 'SpThrSubs_' satstr{sat,1} '_1st4_Qrows.mat']); 
        load([workdir 'SpThrVals_' satstr{sat,1} '_1st4_Qrows.mat']); 
        load([workdir 'SpThrCnts_' satstr{sat,1} '_1st4_Qrows.mat']); 
        
        disp('Processing exclusions.');
        % Origin country exclusions  
        OCEx = {'Cayman Islands'};
        for i=1:size(OCEx,2);
            tmp = find(strcmp(OCEx{1,i},EoraFullLabels(:,1))); 
            tmp1 = ismember(SpThrSubs(:,2),tmp); SpThrSubs = SpThrSubs(~tmp1,:); SpThrVals = SpThrVals(~tmp1,:); SpThrCountries = SpThrCountries(~tmp1,:);
        end; clear OCEx tmp1 tmp2 tmp;
 
        % Destination country exclusions  
        DCEx = {'Cayman Islands'};
        for i=1:size(DCEx,2);
            tmp = find(strcmp(DCEx{1,i},EoraFullLabels(:,1))); 
            tmp2 = ismember(SpThrSubs(:,3),tmp); SpThrSubs = SpThrSubs(~tmp2,:); SpThrVals = SpThrVals(~tmp2,:); SpThrCountries = SpThrCountries(~tmp2,:);
        end; clear DCEx tmp1 tmp2 tmp;
        
        % Producing sector exclusions
        % PSEx = {'Hotels and Restraurants','Construc-tion','Other Constructions','Building and Construction','Civil Engineering','Basic construction','Buildings & constructions','Other Business Activities','Construction ','Construction','Constraction','Other construction','Post and telecommunications','Education, Health and Other Services','Personal Services','Interest groups','Electrical and Machinery','Public Administration','Other business activities','Recycling','Business services','Construction of Communication Facilities','Passenger car trade and repairs, service stations','Life Insurnce Service','Toys and games','Other services','Other Services','other services (private)','Technical services for agriculture, forestry, livestock and fishing','Finacial Intermediation and Business Activities'};
        PSEx = {'Photographic activities and other business activities n.e.c.','Businessactivities','Bussines services except financial services and real estate','Social services','Manufacture of mattresses','Soap','Public adminis-tration and defence; compulsory social security','Social insurance services','Insurance and pension funding, except compulsory social security','Hotels and Restraurants','Machinery & equipment, nec','DOMESTIC SERVICES','Retail Trade','M&E Repair','Activities of membership organisation n.e.c.','Maintenance and Repair','Education, Health and Other Services','Hospital','Electrical and Machinery','Transport Equipment','Interest groups','Other Community Service Activities, Social and Personal Services nec','Post and Telecommunications','Public Administration','Other business activities','Recycling','Business services','Construction of Communication Facilities','Passenger car trade and repairs, service stations','Life Insurnce Service','Toys and games','Other services','Other Services','other services (private)','Technical services for agriculture, forestry, livestock and fishing','Finacial Intermediation and Business Activities'};
        for i=1:size(PSEx,2);
            tmp = find(strcmp(PSEx{1,i},EoraFullLabels(:,2))); 
            tmp1 = ismember(SpThrSubs(:,2),tmp); SpThrSubs = SpThrSubs(~tmp1,:); SpThrVals = SpThrVals(~tmp1,:); SpThrCountries = SpThrCountries(~tmp1,:); 
        end; clear PSEx tmp1 tmp2 tmp;

        % Consuming sector exclusions
        % CSEx = {'Agriculture','Services related to agriculture and forestry','Recycling','Crude oil','Life Insurnce Service','Logging','Agriculture, hunting, forestry and fishing'};
        CSEx = {'Agriculture','Services related to agriculture and forestry','Recycling','Crude oil','Life Insurnce Service','Logging','Agriculture, hunting, forestry and fishing'};
        for i=1:size(CSEx,2);
            tmp = find(strcmp(CSEx{1,i},EoraFullLabels(:,2))); 
            tmp2 = ismember(SpThrSubs(:,3),tmp); SpThrSubs = SpThrSubs(~tmp2,:); SpThrVals = SpThrVals(~tmp2,:); SpThrCountries = SpThrCountries(~tmp2,:);
        end; clear CSEx tmp1 tmp2 tmp;
        
        [tmp,ind] = sort(SpThrVals,'descend'); SpThrSubs = SpThrSubs(ind,:); SpThrCountries = SpThrCountries(ind,:); SpThrVals = tmp;

    end; % skip if results exist clause

    TOPNUMBER = size(SpThrSubs,1);
    RankList = cell(TOPNUMBER,7); 
    RankList(:,1) = SpecName(1,SpThrSubs(1:TOPNUMBER,1)); 
    RankList(:,2) = EoraFullLabels(SpThrSubs(1:TOPNUMBER,2),1); 
    RankList(:,3) = EoraFullLabels(SpThrSubs(1:TOPNUMBER,2),2); 
    RankList(:,4) = EoraFullLabels(SpThrSubs(1:TOPNUMBER,3),1); 
    RankList(:,5) = EoraFullLabels(SpThrSubs(1:TOPNUMBER,3),2); 
    RankList(:,6) = num2cell(SpThrVals(1:TOPNUMBER,1)); 
    RankList(:,7) = countryNames(SpThrCountries(1:TOPNUMBER,1),1);

    save([workdir 'SpThrSubs_' satstr{sat,1} '_final.mat'],'SpThrSubs');
    save([workdir 'SpThrVals_' satstr{sat,1} '_final.mat'],'SpThrVals');
    save([workdir 'SpThrCnts_' satstr{sat,1} '_final.mat'],'SpThrCountries');
    
    save([workdir 'SpThrList_' satstr{sat,1} '_final.mat'],'RankList');
    RankListShort = RankList(1:100000,:);
    save([workdir 'SpThrList_' satstr{sat,1} '_finalShort.mat'],'RankListShort');

    disp(['*** ' satstr{sat,1} ' ***']);
    tmpD = ~strcmp(RankList(:,2),RankList(:,4));
    disp([num2str(sum(cell2mat(RankList(:,6)).*tmpD)) ' international threats; ' num2str(sum(cell2mat(RankList(:,6)))) ' domestic threats.']);
    disp([num2str(sum(cell2mat(RankList(:,6)).*tmpD)./sum(cell2mat(RankList(:,6)))*100) '% international.']);
    disp('.');
        
end; % satellite loop


close all; 
