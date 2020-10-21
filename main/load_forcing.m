%**************************************************************************************************************
% FUNCTION load_forcing.m
% Load forcings preprocessed by "preprocess.m" for the simulation :
% Ecological.mat
% Economical.mat
%**************************************************************************************************************
function forcing = load_forcing(boats,forcing_ecological,forcing_economical,ind_ensemble)

%---------------------------------
% Forcing ecological:
if exist(forcing_ecological,'file')
    load(forcing_ecological);
else
    disp('Hum, double-check the path for ecological forcing:')
    disp(forcing_ecological)
end
forcing.mask=repmat(Ecological.mask,[1 1 size(Ecological.npp,3)]);
forcing.nlat=size(forcing.mask,1);
forcing.nlon=size(forcing.mask,2);
forcing.npp=Ecological.npp;
forcing.npp(find(forcing.mask==1))=NaN;
forcing.npp_ed=Ecological.npp_ed;
forcing.npp_ed(find(forcing.mask==1))=NaN;
forcing.temperature=Ecological.temperature;
forcing.temperature_K=Ecological.temperature+boats.param.conversion.C_2_K;
forcing.surf=Ecological.surface;

%--------------------------------- 
% Forcing economical
if strcmp(boats.param.main.sim_type,'h')
    if exist(forcing_economical,'file')
        load(forcing_economical);
    else
        disp('Hum, double-check the path for economical forcing:')
        disp(forcing_economical)
    end
    load(forcing_economical)
    forcing.price=Economical.price;
    forcing.cost=Economical.cost;
    if (size(Economical.catchability,1) == 1) || (size(Economical.catchability,2) == 1)
        forcing.catchability=Economical.catchability;
    elseif size(Economical.catchability,2) == 5
        forcing.catchability=Economical.catchability(:,ind_ensemble);
    elseif size(Economical.catchability,1) == 5
        forcing.catchability=Economical.catchability(ind_ensemble,1);
    end
    if ~strcmp(boats.param.main.reg_type,'oa')
        if length(size(Economical.effEtarg)) == 5
            forcing.effEtarg=Economical.effEtarg(:,:,:,:,ind_ensemble);             
        elseif length(size(Economical.effEtarg)) == 4
            forcing.effEtarg=Economical.effEtarg(:,:,:,ind_ensemble);             
        end
        forcing.societenf=Economical.societenf;                                
    end
 end


%**************************************************************************************************************
% END FUNCTION

