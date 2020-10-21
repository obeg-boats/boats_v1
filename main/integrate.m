%**************************************************************************************************************
% FUNCTION integrate.m
%**************************************************************************************************************

%-----------------------------------------------------------------------------------------
% Authors
%-----------------------------------------------------------------------------------------

% David A. Carozza  (david.carozza@gmail.com)
% Daniele Bianchi   (dbianchi@atmos.ucla.edu)
% Eric D. Galbraith (eric.d.galbraith@gmail.com)
% Kim J.N. Scherrer (kim.jn.scherrer@gmail.com)

%-----------------------------------------------------------------------------------------
% Introduction
%-----------------------------------------------------------------------------------------

% Bioeconomic Open-Access Trophic Size-Spectrum (BOATS) model based on the
% McKendrick-von Foerster model for a mass-spectrum of fish biomass
% with an open access economic framework to calculate effort and harveSTRU.
% Forced with monthly primary production and temperature data (from model or observations)

% Primary production in units of mmolC m-2 s-1
% Core unit of biomass is mmolC for fish and harvest plots in model
% Convert these to grams of wet biomass (g wetB) using mmolC_2_wetB
% Time in seconds
% Time step is 30 days (1 month)

% 3 fish groups defined by different asymptotic masses
% 1 effort group and 1 selectivity function per fish group
% 2-dimensional version of the code

%-----------------------------------------------------------------------------------------
% Other comments
%-----------------------------------------------------------------------------------------

% We set a lower limit on effort by adding a small constant (epsln = 1e-15)
% to each application of effort when calculating harvest, cost, and effort change.
% This guarantees that as effort approaches zero, the effort change equation approaches 
% the analytical simplification. This prevents dividing by zero.
% We cannot use the analytical simplification of the change in effort directly because we
% limit harvest and we need harvest directly to calculate revenue.yy

% dbianchi 04/16/2016 : output routine changed to allow for general and flexible 
% output "modes". Each mode specifies the variables to be saved, the processing of variables
% (e.g. averages, integrals) and the interval time used for averaging variables (e.g. annual,
% decadal, final timestep, user defined etc.) Additionally, the main calculations were 
% optimized for fast matrix calculations using "bsxfun", and tested with matlab's built-in 
% profiler.

% kscherrer 10/20/2020 : model with dynamic fisheries regulation component 
% as an option to the basic open access model. 

%**************************************************************************************************************
% MAIN CODE
%**************************************************************************************************************
function boats = integrate(boats)

 % start timer
 tic

 %---------------------------------
 % Aliases for variables category for readability
 MAIN=boats.param.main;
 CONV=boats.param.conversion;
 ENVI=boats.param.environment;
 ECOL=boats.param.ecology;
 ECON=boats.param.economy;
 FORC=boats.forcing;
 STRU=boats.structure;
 
 
 %--------------------------------------------------------------------------------------------------------------
 % Defines Ecologic simulations
 %--------------------------------------------------------------------------------------------------------------
 
 %-------------------------------------------------------------------------------------
 % Time integration parameters
 time                   = [MAIN.dtts:MAIN.dtts:MAIN.run_length*CONV.spery]';		   % seconds
 ntime                  = length(time);

 %-----------------------------------------------------------------------------------------
 % Update ecological parameters
 A0                     = ECOL.A00/CONV.spery;                             % growth rate per second (Andersen and Beyer, 2013, p. 4)

 %-----------------------------------------------------------------------------------------
 % Initialize biological arrays
 en_input_P             = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);
 en_input_vb        	= nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);
 en_input               = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);
 gamma                  = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);
 flux_in                = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);
 flux_out               = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);
 flux_fish_growth       = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);
 flux_in_rep            = nan(FORC.nlat,FORC.nlon,ECOL.nfish);
 flux_in_P              = nan(FORC.nlat,FORC.nlon,ECOL.nfish);
 flux_in_num_eggs       = nan(FORC.nlat,FORC.nlon,ECOL.nfish);
 ena_regime             = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass); 
 mortality              = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);

 %-----------------------------------------------------------------------------------------
 % Biomass initial condition
 dfish                  = boats.initial.dfish;
 % Partition of primary production between groups
 part_PP_b              = 1/ECOL.nfish;                                    % partition of primary production at boundary condition (recruitment)
 part_PP_g              = 1/ECOL.nfish;                                    % partition of primary production of growth 
 
 
 %--------------------------------------------------------------------------------------------------------------
 % Defines fileds for Economics simulations
 %--------------------------------------------------------------------------------------------------------------
 if strcmp(MAIN.sim_type,'h')

   %-----------------------------------------------------------------------------------------
   % Initialize economics arrays
   dharvest             = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);
   dfish_temp           = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);      
   gr_prod              = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass); % array to save out gross production
   mort                 = nan(FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass); % array to save out mortality
   effort               = nan(FORC.nlat,FORC.nlon,ECOL.nfish);
   effort_month         = nan(FORC.nlat,FORC.nlon,ECOL.nfish);
   effort_change        = nan(FORC.nlat,FORC.nlon,ECOL.nfish);
   cost                 = nan(FORC.nlat,FORC.nlon,ECOL.nfish);
   revenue              = nan(FORC.nlat,FORC.nlon,ECOL.nfish);
   revenue_memory       = nan(FORC.nlat,FORC.nlon,ECOL.nfish,12);
   cost_memory          = nan(FORC.nlat,FORC.nlon,ECOL.nfish,12);
   effort_memory        = nan(FORC.nlat,FORC.nlon,ECOL.nfish,12);
   regulation_onset     = zeros(FORC.nlat,FORC.nlon,ECOL.nfish); % 0 = no reg target set, 1 = regulation has been initiated
   %-----------------------------------------------------------------------------------------
   % Initialize regulation arrays 
   if strcmp(MAIN.reg_type,'omnT')
       % time series for determining regulation onset
       harvest_sum              = nan(FORC.nlat,FORC.nlon,ECOL.nfish);
       harvest_times            = zeros(FORC.nlat,FORC.nlon,ECOL.nfish,ECON.times_length);
       annual_harvest           = zeros(FORC.nlat,FORC.nlon,ECOL.nfish); % annual total harvest used for economic calculations??
       annual_harvest_max       = zeros(FORC.nlat,FORC.nlon,ECOL.nfish); % annual harvest max used for regulation onset
       regulation_onset         = zeros(FORC.nlat,FORC.nlon,ECOL.nfish); % Needed if only one effort target should be calculated, 0 = no reg target set, 1 = regulation has been initiated,  [* ECON.reg_onset]
       indt_at_reg_onset        = zeros(FORC.nlat,FORC.nlon,ECOL.nfish);
       
   end
   
   %-----------------------------------------------------------------------------------------
   % Initialize some additional output array
   catchability_used 	= nan(1,ntime);
   price_used           = nan(1,ntime);
   cost_effort_used 	= nan(1,ntime);
%    societenf_used       = nan(1,ntime); % ONLY IF Se IS A CONSTANT
%    effEtarg_used        = zeros(FORC.nlat,FORC.nlon,ECOL.nfish,ntime);

   %-----------------------------------------------------------------------------------------
   % Effort initial condition
   effort(:,:,:)         =  boats.initial.effort;
   
 end


 %--------------------------------------------------------------------------------------------------------------
 % Prepare output subroutine
 %--------------------------------------------------------------------------------------------------------------
 outmode = boats.output;
 outmode.noutm = length(outmode.modes);
 % Prepare averaging time bounds
 % This is for cases where time bounds are not specified
 for indm=1:outmode.noutm
    tmode = outmode.modes{indm};
    outmode.(tmode).nvar = length(outmode.(tmode).var_name);
    %----------------------------------------------------------------------
    % Defines the time intervals for averaging
    % here anly the defaults are included
    switch outmode.modes{indm}
    case 'all'
       % NOTE: For output every timestep, no checks/averaging required
       odt = (30*sperd) /dtts; 
       oendt = time(end)/dtts;
       outmode.(tmode).it_bounds = [ [1:odt:oendt]' [odt:odt:oendt]'];
    case 'annual'
       % BOATS assumption is 12 months of 30 days 
       odt = (30*12*CONV.sperd) /MAIN.dtts; 
       oendt = time(end)/MAIN.dtts;
       outmode.(tmode).it_bounds = [ [1:odt:oendt]' [odt:odt:oendt]'];
    case 'decadal'
       odt = (10*30*12*sperd) /dtts; 
       oendt = time(end)/dtts;
       outmode.(tmode).it_bounds = [ [1:odt:oendt]' [odt:odt:oendt]'];
    case 'snap5year'
       % one year every 5 years
       odt  = (5*30*12*sperd) /dtts; 
       odt1 = (1*30*12*sperd) /dtts; 
       oendt = time(end)/dtts;
       vec2 = [odt:odt:oendt]';
       vec1 = vec2-odt1+1;
       outmode.(tmode).it_bounds = [vec1 vec2];
    case 'snap10year'
       % one year every 10 years
       odt  = (10*30*12*sperd) /dtts; 
       odt1 = (1*30*12*sperd) /dtts; 
       oendt = time(end)/dtts;
       vec2 = [odt:odt:oendt]';
       vec1 = vec2-odt1+1;
       outmode.(tmode).it_bounds = [vec1 vec2];
    case 'final'
       % Final year only
       odt = (30*12*sperd) /dtts; 
       oendt = time(end)/dtts;
       outmode.(tmode).it_bounds = [oendt-odt+1 oendt];
    otherwise
       % Uses user-defined time bounds and convert from years
       outmode.(tmode).it_bounds = outmode.(tmode).t_bounds*spery/dtts;
       outmode.(tmode).it_bounds(outmode.(tmode).it_bounds==0) = 1;
    end
    %----------------------------------------------------------------------
    % correct for case when only 1 time integral is given, but it's valid
    outmode.(tmode).t_bounds = time(outmode.(tmode).it_bounds);
    if (size(outmode.(tmode).t_bounds,1)>1)&(size(outmode.(tmode).t_bounds,2)==1); 
       outmode.(tmode).t_bounds = outmode.(tmode).t_bounds';
    end
    outmode.(tmode).time = mean(outmode.(tmode).t_bounds,2);
    outmode.(tmode).year = outmode.(tmode).time/CONV.spery;
    outmode.(tmode).ndt = diff(outmode.(tmode).it_bounds,[],2)+1;
    outmode.(tmode).ntime = length(outmode.(tmode).time);
    % Does a check of variables - if variables have not been defined, removes them from output 
    ivbad = [];   
    for indv=1:outmode.(tmode).nvar
       if ~exist(outmode.(tmode).var_name{indv})
          ivbad = [ivbad indv];
       end
    end
    outmode.(tmode).var_name(ivbad) = [];
    outmode.(tmode).var_type(ivbad) = [];
    outmode.(tmode).var_proc(ivbad) = [];
    outmode.(tmode).var_outn(ivbad) = [];
    outmode.(tmode).nvar = length(outmode.(tmode).var_name);
    if isfield(outmode.(tmode),'var_derv')
       outmode.(tmode).var_derv(ivbad) = [];
    end
 end
 % Removes all invalid output modes:
 % (1) no variables to save
 % (2) no timesteps to save
 imbad = [];
 for indm=1:outmode.noutm
    nvar = outmode.(outmode.modes{indm}).nvar;
    itsz = size(outmode.(outmode.modes{indm}).it_bounds);
    if nvar==0|itsz(2)==1|any(itsz==0)
       imbad = [imbad indm];
    end
 end
 outmode = rmfield(outmode,outmode.modes(imbad));
 outmode.modes(imbad) = [];
 outmode.noutm = length(outmode.modes);

 
 %--------------------------------------------------------------------------------------------------------------
 % Initialize output arrays
 %--------------------------------------------------------------------------------------------------------------
 % This require hardwiring the dimensions of output arrays based on the
 % variable original size/type (4D,3D, etc) and processing (2di, si, gi)
 disp(['Initializing output...']);
 for indm=1:outmode.noutm
    ontime = outmode.(outmode.modes{indm}).ntime;
    for indv=1:outmode.(outmode.modes{indm}).nvar
      % Defines output size according to variable type
       switch outmode.(outmode.modes{indm}).var_type{indv}
       case '4D'
          % Further defines output size according to variable processing
          switch outmode.(outmode.modes{indm}).var_proc{indv}
          case 'none'
          % Saves the entire variable without any processing
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,ECOL.nfish,ECOL.nfmass);
          case 'si'
          % Saves the variable after integral over size dimensions
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,ECOL.nfish);
          case '2di_si'
          % Saves the variable after integral over size dimensions and over 2D domain
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,ECOL.nfish);
          case 'LMEi_si'
          % Saves the variable after integral over size dimensions and over 2D domain only across LME
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,ECOL.nfish);
          case '2di_si_gi'
          % Saves the variable summed over the fish groups after integral over size dimensions and over 2D
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,1);
          otherwise
             error(['processing case ' outmode.(outmode.modes{indm}).var_proc{indv} ' not specified']);
          end
       case '3D'
          % Further defines output size according to variable processing
          switch outmode.(outmode.modes{indm}).var_proc{indv}
          case 'none'
          % Saves the entire variable without any processing
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,FORC.nlat,FORC.nlon,ECOL.nfish);
          case '2di'
          % Saves the variable after integral over 2D domain
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,ECOL.nfish);
          case 'LMEi'
          % Saves the variable after integral over 2D domain only across LME
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,ECOL.nfish);
          case '2di_gi'
          % Saves the variable summed over groups after integral over 2D domain 
             outmode.(outmode.modes{indm}).(outmode.(outmode.modes{indm}).var_outn{indv}) = ...
                                           zeros(ontime,1);
          otherwise
             error(['processing case ' outmode.(outmode.modes{indm}).var_proc{indv} ' not specified']);
          end
       case 'DER'
          % Initialization generally not be required for derived variables
       otherwise
          error(['processing case ' outmode.(outmode.modes{indm}).var_type{indv} ' not specified']);
       end
    end
 end

 
%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------
% MAIN LOOP
%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------

for indt = 1:ntime
    
    % local month
    local_month = ceil((mod(time(indt),CONV.spery)/CONV.spery)*12);
    if local_month==0;local_month=12;end
    disp(['indt : ' num2str(indt) ' / ' num2str(ntime) ' - local month : ' num2str(local_month)]);
    
    %---------------------------------------------------------------------------------------
    % Physical and ecological forcings (NPP in mmolC m-2 s-1)
    %---------------------------------------------------------------------------------------
    if size(FORC.npp,3) >= ntime % For time-varying Ecological forcing  || 2019-12-11: KIM CHANGED CONDITION FROM == to >= TO BE ABLE TO RUN SIMULATIONS ONLY FOR PARTS OF THE FORCINGS TIMES
        npp         = squeeze(FORC.npp(:,:,indt));
        npp_ed      = squeeze(FORC.npp_ed(:,:,indt));
        temp_phyto  = squeeze(FORC.temperature(:,:,indt));
        temp_fish_A = squeeze(FORC.temperature_K(:,:,indt));
    elseif size(FORC.npp,3) == 12 % For monthly climatology forcing
        npp         = squeeze(FORC.npp(:,:,local_month));
        npp_ed      = squeeze(FORC.npp_ed(:,:,local_month));
        temp_phyto  = squeeze(FORC.temperature(:,:,local_month));
        temp_fish_A = squeeze(FORC.temperature_K(:,:,local_month));
    end
    temp_fish_m = temp_fish_A;
    
    %-----------------------------------------------------------------------------------------------------------
    % Large fraction of phytoplankton and representative phytoplankton mass (Dunne)
    %-----------------------------------------------------------------------------------------------------------
    s_over_p = ( -1.0 + ( 1.0 + 4.0 .* npp_ed ./ (exp(ENVI.kappa_eppley.*temp_phyto) .* ...
        ENVI.Prod_star) ).^0.5) .* 0.5;
    frac_lg_du = s_over_p ./ (1.0 + s_over_p);                               % large fraction of PP as in Dunne et al. (2005)
    mphyto = (ENVI.mc_phy_l.^frac_lg_du) .* (ENVI.mc_phy_s.^(1.0 - frac_lg_du));
    
    %-----------------------------------------------------------------------------------------------------------
    % Growth rate
    %-----------------------------------------------------------------------------------------------------------
    % growth rate = production distribution * mass / biomass distribution
    % multiply by growth partition function (part_PP_g)
    %-------------------------------------------------------------------------------------
    % % Original code - changed to optimize calculation
    % en_input_P = repmat(npp./mphyto,[1 1 nfish nfmass]) .* ...
    %              fmass_4d ./ squeeze(dfish + epsln) * part_PP_g;
    % Optimized code by using "bsxfun" instead of repmat
    
    % No iron-lim HNLC (ECOL.tro_sca = scalar)
    %   en_input_P = bsxfun(@times,npp./mphyto, ...
    %                (bsxfun(@rdivide,permute(STRU.fmass_2d.^(ECOL.tro_sca-1),[3 4 1 2]),mphyto.^(ECOL.tro_sca-1)) .* ...
    %                STRU.fmass_4d ./ squeeze(dfish + CONV.epsln) * part_PP_g));
    % With iron-lim HNLC (ECOL.tro_sca = map)
    en_input_P = bsxfun(@times,npp./mphyto, ...
        (bsxfun(@rdivide,permute(repmat(STRU.fmass_2d,[1 1 180 360]),[3 4 1 2]).^ ...
        (repmat(ECOL.tro_sca,[1 1 ECOL.nfish ECOL.nfmass])-1),mphyto.^(ECOL.tro_sca-1)) .* ...
        STRU.fmass_4d ./ squeeze(dfish + CONV.epsln) * part_PP_g));
    
    %---------------------------------------------------------------------------------------
    % Based on allometric scaling (von Bertalanffy)
    % calculate temperature dependencies, then growth rate, the activity loss constant (ka)
    temp_dep_A  = repmat(exp( (-ENVI.E_activation_A/ENVI.k_Boltzmann) .* (1./temp_fish_A - 1/ENVI.temp_ref_A)),[1 1 ECOL.nfish ECOL.nfmass]);
    temp_dep_m  = repmat(exp( (-ENVI.E_activation_m/ENVI.k_Boltzmann) .* (1./temp_fish_m - 1/ENVI.temp_ref_A)),[1 1 ECOL.nfish ECOL.nfmass]);
    A 	      = A0 .* temp_dep_A;
    ka 	      = A * ECOL.eff_a .* STRU.minf_4d_p_bm1;                      %  A * eff_a .* minf_4d.^(b_allo-1) (s-1)
    en_input_vb = A .* STRU.fmass_4d_p_b - ka .* STRU.fmass_4d;              % A .* fmass_4d.^b_allo - ka .* fmass_4d;
    
    %---------------------------------------------------------------------------------------
    % Input energy (energy available to growth and reproduction)
    % minimum of en_input_P and en_input_vb
    en_input                     = min(en_input_P,en_input_vb);
    en_input(STRU.mask_notexist_4d)   = NaN;
    
    %---------------------------------------------------------------------------------------
    % Somatic growth rate (growth rate of body tissues)
    gamma                        = (1 - STRU.rep_alloc_frac).*en_input;
    
    %---------------------------------------------------------------------------------------
    % Flux out of a mass class
    flux_out   = gamma .* squeeze(dfish+CONV.epsln) ./ STRU.delfm_4d;
    
    %---------------------------------------------------------------------------------------
    % Biomass growth
    % increase in biomass in growth to next mass class
    % arises in conversion from abundance-based MFVF to biomass-based
    flux_fish_growth   = gamma .* squeeze(dfish+CONV.epsln) ./ STRU.fmass_4d;
    
    %---------------------------------------------------------------------------------------
    % boundary condition (flux in to first mass class)
    %---------------------------------------------------------------------------------------
    % Boundary condition based on primary production
    % multiply by boundary condition partition function (part_PP_b)
    % No iron-lim HNLC (ECOL.tro_sca = scalar)
    %   flux_in_P = part_PP_b * (repmat(npp./mphyto,[1 1 ECOL.nfish])) .* (STRU.fmass_bc./repmat(mphyto,[1 1 ECOL.nfish])).^ ...
    %               (ECOL.tro_sca-1) * STRU.fmass_bc / STRU.delfm_4d(1,1,1,1);
    % With iron-lim HNLC (ECOL.tro_sca = map)
    flux_in_P = part_PP_b * (repmat(npp./mphyto,[1 1 ECOL.nfish])) .* (STRU.fmass_bc./repmat(mphyto,[1 1 ECOL.nfish])).^ ...
        (repmat(ECOL.tro_sca,[1 1 ECOL.nfish])-1) * STRU.fmass_bc / STRU.delfm_4d(1,1,1,1);
    
    %---------------------------------------------------------------------------------------
    % Flux in of number of eggs produced
    flux_in_num_eggs = (ECOL.frac_fem/ECOL.m_egg) .* ...
        nansum( STRU.rep_alloc_frac .* en_input .* squeeze(dfish) .* STRU.delfm_4d_over_fmass_4d,4) / STRU.delfm_4d(1,1,1,1);
    
    %---------------------------------------------------------------------------------------
    % Boundary condition based on recruitment (production and survival of eggs)
    % If the flux of number of eggs is less than 0.001 eggs m-2 y-1
    % then the flux due to the production and survival of eggs is zero
    mask_eggs_low = (flux_in_num_eggs < 0.001/CONV.spery);
    flux_in_rep = (STRU.fmass_bc*ECOL.egg_surv) .* flux_in_num_eggs;
    flux_in_rep(mask_eggs_low) = 0;
    
    %---------------------------------------------------------------------------------------
    % Boundary condition (Beverton-Holt form)
    flux_in(:,:,:,1) = flux_in_P .* ((flux_in_rep + CONV.epsln) ./ (flux_in_P + flux_in_rep + CONV.epsln));
    
    %---------------------------------------------------------------------------------------
    % Flux in of other mass classes
    flux_in(:,:,:,2:end) = gamma(:,:,:,1:end-1) .* squeeze(dfish(:,:,:,1:end-1)) ./ STRU.delfm_2end_4d;
    
    %---------------------------------------------------------------------------------------
    % Input energy (available energy to growth and reproduction regime
    % 1 if en_input_P less; then en_input_P determines input energy
    ena_regime = (en_input_P < en_input_vb);
    
    %---------------------------------------------------------------------------------------
    % Mortality
    %---------------------------------------------------------------------------------------
    % Mass-specific mortality rate
    % Calculate associated growth rate with mortality temperature dependence temp_dep_m
    % calculate mortality rate mortality0
    % mortality00 is the exp(zeta_1) term from the model description
    A           = A0*temp_dep_m;
    mortality0  = (exp(ECOL.zeta1)/3)*A;
    
    %---------------------------------------------------------------------------------------
    % Mortality rate
    % Charnov et al. (2013)
    % minf_4d_p_hplusbm1 is minf_4d.^(h_allo + b_allo - 1)
    mortality = mortality0 .* STRU.minf_4d_p_hplusbm1 .* STRU.fmass_4d_p_mh .* squeeze(dfish);
    
    %---------------------------------------------------------------------------------------
    % Integrate fish and effort
    %---------------------------------------------------------------------------------------
    if strcmp(MAIN.sim_type,'nh') % No economic harvesting
        
        %-------------------------------------------------------------------------------------
        % Integrate dfish
        %-------------------------------------------------------------------------------------
        dfish  = squeeze(dfish) + ( flux_in - flux_out + flux_fish_growth - mortality ) * MAIN.dtts;
        mask_dfish_neg = (squeeze(dfish) < 0);
        dfish(mask_dfish_neg)  = 0;
        
    else % Economic harvesting
        
        %-------------------------------------------------------------------------------------
        % Update fish to calculate harvest
        %-------------------------------------------------------------------------------------
        dfish_temp  = squeeze(dfish) + ( flux_in - flux_out + flux_fish_growth - mortality ) * MAIN.dtts;
        mask_dfish_temp_neg = (squeeze(dfish_temp) < 0);
        dfish_temp(mask_dfish_temp_neg)  = 0;
        
        %Save out eggs, mortality, fish production
        gr_prod = flux_in - flux_out + flux_fish_growth;
        mort = mortality;
        gr_prod(mask_dfish_temp_neg) = mortality(mask_dfish_temp_neg) - squeeze(dfish(mask_dfish_temp_neg)) / MAIN.dtts;
        % Determine values in cells with negative dfish;
        % gr_prod = -dfish(lasttimestep)/dtts + mortality
        
        % flux_in - flux_out + flux_fish_growth = fish production
        
        
        %-------------------------------------------------------------------------------------
        % Catchability
        %-------------------------------------------------------------------------------------
        if (time(indt) < ECON.harvest_start*CONV.spery)
            % Catchability zero before start year (adjusted for dharvest calculation)
            qcatch = 0 * ones(1,ECOL.nfish);
        else
            % Set catchability forcing scenario (adjusted for dharvest calculation)
            catchability_used(1,indt) = FORC.catchability(indt);
            qcatch = FORC.catchability(indt) * ones(1,ECOL.nfish);
        end
        
        %-------------------------------------------------------------------------------------
        % dharvest [nlat,nlon,nfish,nfmass]
        %-------------------------------------------------------------------------------------
        % qcatch * effort * selectivity * dfish_temp
        % Set upper limit (min) so that no more than the fish present can be caught (dfish_temp/dtts)
        %-------------------------------------------------------------------------------------
        % Original code - changed to optimize calculation
        % dharvest = min(squeeze(dfish_temp)/dtts, permute(repmat(qcatch(:),[1 nlat nlon nfmass]),[2 3 1 4]) .* ...
        %      repmat(squeeze(effort+epsln),[1 1 1 nfmass]) .* selectivity_4d .* squeeze(dfish_temp));
        % Optimized code by using "bsxfun" instead of repmat
        dharvest = min(squeeze(dfish_temp)/MAIN.dtts, ...
            bsxfun(@times,bsxfun(@times,permute(qcatch(:),[2 3 1]),squeeze(effort+CONV.epsln)), ... 
            STRU.selectivity_4d .* squeeze(dfish_temp)));
        
        mask_dharvest_neg = (squeeze(dharvest) < 0);
        dharvest(mask_dharvest_neg)  = 0;
        
        %----------------------------------------------------------------------------------
        % Price forcing
        % input price ($ g-1), multiply by mmolC_2_wetB to get ($ mmolC-1)
        %----------------------------------------------------------------------------------
        % Set price forcing scenario (adjusted for dharvest calculation)
        price_used(1,indt) = FORC.price(indt);
        price       = FORC.price(indt) * ones(ECOL.nfish,ECOL.nfmass);
        
        %-------------------------------------------------------------------------------------
        % Cost per unit effort forcing
        % input cost per unit effort ($ W-1 s-1)
        %-------------------------------------------------------------------------------------
        % Set cost forcing scenario (adjusted for dharvest calculation)
        cost_effort_used(1,indt) = FORC.cost(indt);
        cost_effort = FORC.cost(indt) * ones(1,ECOL.nfish);
        
        %-------------------------------------------------------------------------------------
        % revenue [nlat,nlon,nfish]
        %-------------------------------------------------------------------------------------
        % sum over mass (a * dharvest * delfm)
        %-------------------------------------------------------------------------------------
        % Original code - changed to optimize calculation
        % revenue = nansum( permute(repmat(price,[1 1 nlat nlon]),[3 4 1 2]) .* squeeze(dharvest) .* delfm_4d, 4);
        % Optimized code by using "bsxfun" instead of repmat
        revenue = nansum(bsxfun(@times,permute(price,[3 4 1 2]),squeeze(dharvest) .* STRU.delfm_4d),4);
        
        %-------------------------------------------------------------------------------------
        % revenue_memory [nlat,nlon,nfish,12] % SAVES REVENUE FR. LAST 12 MONTHS
        %-------------------------------------------------------------------------------------
        % Save revenue time series to use for annual effort update
        if indt == 1  % set first value of harvest_sum as the first point in time series
            revenue_memory = revenue;
        elseif indt <= 12 % add present harvest_sum value as the following point until reaching desired timeseries lenght (times_length)
            revenue_memory = cat(4, revenue_memory, revenue);  % add new point to the 4th dimension (i.e. the time series) in harvest_data
        else        % when timeseries length is times_length, cut off the oldest data point and add present harvest_sum to the end of time serie
            revenue_memory(:,:,:,1:11) = revenue_memory(:,:,:,2:12);
            revenue_memory(:,:,:,12) = revenue;
        end
        
        %-------------------------------------------------------------------------------------
        % cost [nlat,nlon,nfish]
        %-------------------------------------------------------------------------------------
        % cost_effort * effort
        %-------------------------------------------------------------------------------------
        % Original code - changed to optimize calculation
        % cost =  permute(repmat(cost_effort(:),[1 nlat nlon]),[2 3 1]) .* squeeze(effort + epsln);
        % Optimized code by using "bsxfun" instead of repmat
        cost =  bsxfun(@times,permute(cost_effort(:),[2 3 1]),squeeze(effort + CONV.epsln)); 
        
        %-------------------------------------------------------------------------------------
        % cost_memory [nlat,nlon,nfish,12] % SAVES COST FR. LAST 12 MONTHS
        %-------------------------------------------------------------------------------------
        % Save cost time series to use for annual effort update
        if indt == 1  % set first value of harvest_sum as the first point in time series
            cost_memory = cost;
        elseif indt <= 12 % add present harvest_sum value as the following point until reaching desired timeseries lenght (times_length)
            cost_memory = cat(4, cost_memory, cost);  % add new point to the 4th dimension (i.e. the time series) in harvest_data
        else        % when timeseries length is 12 months long, cut off the oldest data point and add present harvest_sum to the end of time serie
            cost_memory(:,:,:,1:11) = cost_memory(:,:,:,2:12);
            cost_memory(:,:,:,12) = cost;
        end
        
        if strcmp(MAIN.reg_type,'oa') % no regulation
            
            %-------------------------------------------------------------------------------------
            % effort_change [nlat,nlon,nfish] 
            % (monthly)
            %-------------------------------------------------------------------------------------
            effort_change = ECON.k_e * (revenue - cost) ./ (effort + CONV.epsln);
            
            %-------------------------------------------------------------------------------------
            % integrate effort [nlat,nlon,nfish]
            %-------------------------------------------------------------------------------------
            % Monthly integration of effort
            effort = squeeze(effort) + effort_change * MAIN.dtts; 
            mask_effort_neg = (squeeze(effort) < 0);
            effort(mask_effort_neg)  = 0;
            
        elseif strcmp(MAIN.reg_type,'omnT') % regulation
            
            %----------------------------------------------------------------------------------
            % Societal enforcement forcing [nlat,nlon,nfish]
            %----------------------------------------------------------------------------------
            % Set societal enforcement forcing scenario
            if length(size(FORC.societenf)) == 2
                % Se CONSTANT OVER TIME
                societenf = cat(3,FORC.societenf,FORC.societenf,FORC.societenf); % create identical maps for all three groups
            elseif length(size(FORC.societenf)) == 3
                % Se HAS TIME DIMENSION
                societenf = cat(3,squeeze(FORC.societenf(:,:,indt)),squeeze(FORC.societenf(:,:,indt)),squeeze(FORC.societenf(:,:,indt))); % create identical maps for all three groups
            end
            
            %----------------------------------------------------------------------------------
            % Effective effort target forcing [nlat,nlon,nfish]
            %----------------------------------------------------------------------------------
            % Set effective effort target forcing scenario
            if length(size(FORC.effEtarg)) == 3 % effective effort target has no time dimension [lat lon group]
                effEtarg = FORC.effEtarg; % If effEtarg = constant over time
                
            elseif length(size(FORC.effEtarg)) == 4  % indt >= 12 %
                effEtarg = squeeze(FORC.effEtarg(:,:,indt/12,:));
            end
            
            %-------------------------------------------------------------------------------------
            % Regulated and unregulated cells [nlat,nlon,nfish]
            %-------------------------------------------------------------------------------------
            % Calculate effort evolution in cells where a regulation target has been
            % estimated because of decreasing catch
            mask_reg            = find(regulation_onset == 1); % cells where regulation has started
            mask_noreg          = find(regulation_onset == 0);
            % mask_reg           = find((regulation_onset == 1) & (effEtarg * ECON.precaut / catchability_y_avg) < effort ); % cells where regulation has started AND current nominal effort is higher than target nominal effort
            
            %-------------------------------------------------------------------------------------
            % Calculate effort_change every month
            %-------------------------------------------------------------------------------------
            if strcmp(ECON.reg_timestep,'monthly')
                
                %-------------------------------------------------------------------------------------
                % Calculate regulated effort_change
                %-------------------------------------------------------------------------------------
                effort_change(mask_reg) = ECON.k_e * (revenue(mask_reg) - cost(mask_reg)) ./ (effort(mask_reg) + CONV.epsln) ...
                    .* (ones(length(mask_reg),1) - societenf(mask_reg)) + ECON.k_s .* societenf(mask_reg) .* (effEtarg(mask_reg) * ECON.precaut ...
                    / FORC.catchability(indt) - effort(mask_reg)); 
                
                %-------------------------------------------------------------------------------------
                % Calculate unregulated effort_change [nlat,nlon,nfish]
                %-------------------------------------------------------------------------------------
                % Calculate effort evolution in cells where no regulation target has
                % been estimated, i.e. open access effort
                effort_change(mask_noreg)   = ECON.k_e .* (revenue(mask_noreg) - cost(mask_noreg)) ./ (effort(mask_noreg) + CONV.epsln);
                
                %-------------------------------------------------------------------------------------
                % integrate effort [nlat,nlon,nfish]
                %-------------------------------------------------------------------------------------
                % Monthly integration of effort
                effort = squeeze(effort) + effort_change * MAIN.dtts;
                mask_effort_neg = (squeeze(effort) < 0);
                effort(mask_effort_neg)  = 0;
                
                
                %-------------------------------------------------------------------------------------
                % Calculate effort_change once a year
                %-------------------------------------------------------------------------------------
            elseif (mod(indt,12) == 1) && (indt > 12) % reg_timestep == 'yearly', mod(indt,12) == 1 corresponds to January
                
                
                %-------------------------------------------------------------------------------------
                % revenue_y, cost_y, effort_y [nlat,nlon,nfish]
                %-------------------------------------------------------------------------------------
                % Annual sum of rates of revenue, profit, effort and catchability for calculation
                % of effort_change (based on profit-effort ratio and effort-effort
                % target difference).
                
                revenue_y_avg       = mean(revenue_memory,4);
                cost_y_avg          = mean(cost_memory,4);
                effort_y_avg        = mean(effort_memory,4);
                catchability_y_avg  = mean(catchability_used(1,indt-12:indt-1));
                
                %-------------------------------------------------------------------------------------
                % Calculate regulated effort_change
                %-------------------------------------------------------------------------------------
                effort_change(mask_reg) = ECON.k_e * (revenue_y_avg(mask_reg) - cost_y_avg(mask_reg)) ./ (effort_y_avg(mask_reg) + CONV.epsln) ...
                    .* (ones(length(mask_reg),1) - societenf(mask_reg)) + ECON.k_s .* societenf(mask_reg) .* (effEtarg(mask_reg) * ECON.precaut ...
                    / catchability_y_avg - effort_y_avg(mask_reg));
                
                %-------------------------------------------------------------------------------------
                % Calculate unregulated effort_change [nlat,nlon,nfish]
                %-------------------------------------------------------------------------------------
                % Calculate effort evolution in cells where no regulation target has
                % been estimated, i.e. open access effort
                effort_change(mask_noreg)   = ECON.k_e .* (revenue_y_avg(mask_noreg) - cost_y_avg(mask_noreg)) ./ (effort_y_avg(mask_noreg) + CONV.epsln);
                
                %-------------------------------------------------------------------------------------
                % integrate effort [nlat,nlon,nfish]
                %-------------------------------------------------------------------------------------
                % Annual integration of effort
                effort = squeeze(effort) + effort_change * CONV.spery; 
                mask_effort_neg = (squeeze(effort) < 0);
                effort(mask_effort_neg)  = 0;
            end
            
            %-------------------------------------------------------------------------------------
            % effort_memory [nlat,nlon,nfish,12]
            %-------------------------------------------------------------------------------------
            % Save last 12 months of effort to use in calculation of effort_change
            
            if indt == 1  % set first value of harvest_sum as the first point in time series
                effort_memory = effort;
            elseif indt <= 12 % add present harvest_sum value as the following point until reaching desired timeseries lenght (times_length)
                effort_memory = cat(4, effort_memory, effort);  % add new point to the 4th dimension (i.e. the time series) in harvest_data
            else        % when timeseries length is times_length, cut off the oldest data point and add present harvest_sum to the end of time serie
                effort_memory(:,:,:,1:11) = effort_memory(:,:,:,2:12);
                effort_memory(:,:,:,12) = effort;
            end
            
        end
        
        %-------------------------------------------------------------------------------------
        % integrate dfish [nlat,nlon,nfish,nfmass]
        %-------------------------------------------------------------------------------------
        dfish  = squeeze(dfish_temp) - squeeze(dharvest) * MAIN.dtts;
        mask_dfish_neg = (squeeze(dfish) < 0);
        dfish(mask_dfish_neg)  = 0;
        
        
        
        %-------------------------------------------------------------------------------------
        % DETERMINE IF REGULATION SHOULD BEGIN: regulation_onset
        %-------------------------------------------------------------------------------------
        if ~strcmp(MAIN.reg_type,'oa')                                         % if not an OA simulation
            
            %-------------------------------------------------------------------------------------
            % harvest_sum [nlat,nlon,nfish]
            %-------------------------------------------------------------------------------------
            % sum harvest over mass classes (dharvest * delfm) to be used as
            % indicators of need for regulation
            %-------------------------------------------------------------------------------------
            harvest_sum(:,:,:) = nansum(dharvest .* STRU.delfm_4d,4);
            
            %-------------------------------------------------------------------------------------
            % harvest_times [nlat,nlon,nfish,times_length]
            %-------------------------------------------------------------------------------------
            % save "timeseries" of harvest for determination of regulation onset
            %-------------------------------------------------------------------------------------
            if indt == 1  % set first value of harvest_sum as the first point in time series
                harvest_times = harvest_sum;    
            elseif indt <= ECON.times_length % add present harvest_sum value as the following point until reaching desired timeseries lenght (times_length)
                harvest_times = cat(4, harvest_times, harvest_sum);  % add new point to 4th dimension (i.e. the time series) in harvest_data
            else        % when timeseries length is times_length, cut off the oldest data point and add present harvest_sum to the end of time serie
                harvest_times(:,:,:,1:ECON.times_length-1) = harvest_times(:,:,:,2:ECON.times_length);
                harvest_times(:,:,:,ECON.times_length) = harvest_sum;
            end
            
            %-------------------------------------------------------------------------------------
            % annual_harvest_max [nlat,nlon,nfish]
            %-------------------------------------------------------------------------------------
            % save the highest recorded annual harvest so far/of all time steps
            %-------------------------------------------------------------------------------------
            if indt == 1
                annual_harvest_max = harvest_sum; % default initial harvest max is harvest_sum in January 
            elseif mod(indt,12) == 1 % each january, sum the 12 previous months in harvest_data to get annual harvest
                annual_harvest = nansum(harvest_times(:,:,:,end-12:end-1),4); % save harvest from previous 12 months using the last 12 months in harvest_data
                annual_max_ind = find(annual_harvest > annual_harvest_max); % find all cells where calculated annual_harvest is higher than annual_harvest_max previously calculated
                annual_harvest_max(annual_max_ind) = annual_harvest(annual_max_ind); % update new highest annual harvest so far
            end
            
            %-------------------------------------------------------------------------------------
            % effort_times [nlat,nlon,nfish,times_length]
            %-------------------------------------------------------------------------------------
            % save "timeseries" of effort for 'optS' and 'stockT' regulation
            %-------------------------------------------------------------------------------------
            if indt == 1
                effort_times           = effort; 
            elseif indt <= ECON.times_length
                effort_times           = cat(4, effort_times, effort); 
            else
                effort_times(:,:,:,1:ECON.times_length-1)  = effort_times(:,:,:,2:ECON.times_length);
                effort_times(:,:,:,ECON.times_length)      = effort; 
            end
            
            %-------------------------------------------------------------------------------------
            % DETERMINE REGULATION ONSET (AND TARGET DEPENDING ON reg_type)
            %-------------------------------------------------------------------------------------
            for lat = 1:FORC.nlat
                for lon = 1:FORC.nlon
                    for group = 1:ECOL.nfish
                        if (mod(indt,12) == 1) && (indt ~= 1) && (annual_harvest(lat,lon,group) < ECON.reg_threshold * ...
                                annual_harvest_max(lat,lon,group)) && (regulation_onset(lat,lon,group) == 0) ...
                                && (annual_harvest(lat,lon,group) > 1e-25 ) % Every January, determine regulation onset and target, last condition added to avoid too-early regulation onset
                            if strcmp(ECON.reg_timestep,'monthly')
                                % Set regulation onset to 1 and start the adaptive regulation in next step
                                if ((effEtarg(lat,lon,group) * ECON.precaut / FORC.catchability(indt)) < effort(lat,lon,group))
                                    regulation_onset(lat,lon,group) = 1;
                                    % save the 'indt' for time of regulation onset
                                    indt_at_reg_onset(lat,lon,group) = indt;
                                end
                            elseif strcmp(ECON.reg_timestep,'yearly')
                                if ((effEtarg(lat,lon,group) * ECON.precaut / catchability_y_avg) < effort(lat,lon,group))
                                    regulation_onset(lat,lon,group) = 1;
                                    % save the 'indt' for time of regulation onset
                                    indt_at_reg_onset(lat,lon,group) = indt;
                                end
                            end
                        end
                    end
                end
            end
        end
    end % end integrate
    
    
%-----------------------------------------------------------------------------------------  
% Output : calculate and save output variables
%-----------------------------------------------------------------------------------------
% Fills in output variables inside previously defined output "modes" structures. 
% Handles the time averaging dynamically by dividing the temporal output by the number 
% of timesteps and summing. Sadly, this part requires a few "eval statements"
 for indm=1:outmode.noutm
    % Decides whether the current timestep should be saved out 
    % finds the index corresponding to the current time averaging interval - if it exists
    prodt = prod(outmode.(outmode.modes{indm}).it_bounds - indt,2);
    iit = find(prodt<=0,1,'first');
    if ~isempty(iit);
       % Processes the current output timestep
       for indv=1:outmode.(outmode.modes{indm}).nvar
          tvar_outn = outmode.(outmode.modes{indm}).var_outn; 
          % Defines output size according to variable type
          switch outmode.(outmode.modes{indm}).var_type{indv}
          case '4D'
             tmpvar = eval(outmode.(outmode.modes{indm}).var_name{indv});
             switch outmode.(outmode.modes{indm}).var_proc{indv}
             case 'none'
             % Saves the entire variable without any processing
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:,:))  + ...
                                              tmpvar/outmode.(outmode.modes{indm}).ndt(iit) + ...
                                              STRU.mask_land_g_s_nan;
             case 'si'
             % Saves the variable after integral over size dimensions
                tmpvar1 = squeeze(nansum( tmpvar .* STRU.delfm_4d,4)) * CONV.mmolC_2_wetB + STRU.mask_land_g_nan;
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:))  + ...
                                              tmpvar1/outmode.(outmode.modes{indm}).ndt(iit);
             case '2di_si'
             % Saves the variable after integral over size dimensions and over 2D domain
                tmpvar1 = squeeze(nansum(tmpvar .* STRU.delfm_4d,4)) * CONV.mmolC_2_wetB;
                tmpvar2 = nan([1 3]);
                for indgg=1:3
                   tmpvar2(1,indgg) = integrate_2d(squeeze(tmpvar1(:,:,indgg)),boats.forcing.surf);
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar2/outmode.(outmode.modes{indm}).ndt(iit);
             case 'LMEi_si'
             % Saves the variable after integral over size dimensions and over 2D domain only across LME
                tmpvar1 = squeeze(nansum(tmpvar .* STRU.delfm_4d,4)) * CONV.mmolC_2_wetB;
                tmpvar2 = nan([1 3]);
                for indgg=1:3
                   tmpvar2(1,indgg) = integrate_2d(mask_LME_nan+squeeze(tmpvar1(:,:,indgg)),boats.surf);
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar2/outmode.(outmode.modes{indm}).ndt(iit);
             case '2di_si_gi'
             % Saves the variable summed over the fish groups after integral over size dimensions and over 2D
                tmpvar1 = squeeze(nansum(tmpvar .* delfm_4d,4)) * mmolC_2_wetB;
                tmpvar2 = nan([1 3]);
                for indgg=1:3
                   tmpvar2(1,indgg) = integrate_2d(squeeze(tmpvar1(:,:,indgg)),boats.surf);
                end
                tmpvar3 = nansum(tmpvar2,2);
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit))  + ...
                                              tmpvar3/outmode.(outmode.modes{indm}).ndt(iit);
             otherwise
                error(['processing case ' outmode.(outmode.modes{indm}).var_proc{indv} ' not specified']);
             end
          case '3D'
             tmpvar = eval(outmode.(outmode.modes{indm}).var_name{indv});
             % Further defines output size according to variable processing
             switch outmode.(outmode.modes{indm}).var_proc{indv}
             case 'none'
             % Saves the entire variable without any processing
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:,:,:))  + ...
                                              tmpvar/outmode.(outmode.modes{indm}).ndt(iit) + ...
                                              STRU.mask_land_g_nan;
             case '2di'
             % Saves the variable after integral over 2D domain
                tmpvar1 = nan([1 3]);
                for indgg=1:3
                   tmpvar1(1,indgg) = integrate_2d(squeeze(tmpvar(:,:,indgg)),boats.forcing.surf);
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar1/outmode.(outmode.modes{indm}).ndt(iit);
             case 'LMEi'
             % Saves the variable after integral over 2D domain only across LME
                tmpvar1 = nan([1 3]);
                for indgg=1:3
                   tmpvar1(1,indgg) = integrate_2d(mask_LME_nan+squeeze(tmpvar(:,:,indgg)),boats.surf);
                end
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar1/outmode.(outmode.modes{indm}).ndt(iit);
             case '2di_gi'
             % Saves the variable summed over groups after integral over 2D domain
                tmpvar1 = nan([1 3]);
                for indgg=1:3
                   tmpvar1(1,indgg) = integrate_2d(squeeze(tmpvar(:,:,indgg)),boats.surf);
                end
                tmpvar2 = nansum(tmpvar1,2);
                outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:)  = ...
                squeeze(outmode.(outmode.modes{indm}).(tvar_outn{indv})(iit,:))  + ...
                                              tmpvar2/outmode.(outmode.modes{indm}).ndt(iit);
             otherwise
                error(['processing case ' outmode.(outmode.modes{indm}).var_proc{indv} ' not specified']);
             end
          case 'DER'
             % Calculations for derived variables may not be required
          otherwise
             error(['processing case ' outmode.(outmode.modes{indm}).var_type{indv} ' not specified']);
          end
       end
    end
 end
end % for indt = 1:ntime

%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------
% END OF MAIN LOOP
%--------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------


%--------------------------------------------------------------------------------------------------------------
% Save last timestep of dfish and effort for restart
%--------------------------------------------------------------------------------------------------------------
 if MAIN.save_restart==1
   boats.dfish      = dfish;
   if strcmp(MAIN.sim_type,'h')
     boats.effort   = effort; 
   end
 end

%-----------------------------------------------------------------------------------------
% Output : additional derived variables needed for output - 
% these are calculated outside the time-stepping loop
%-----------------------------------------------------------------------------------------
 for indm=1:outmode.noutm
    % Derived variables for output
    % Generally avoid doing calculation and integrals within 
    % the main  time loop if it can be avoided
    ider = find(strcmp(outmode.(outmode.modes{indm}).var_type,'DER'));
    for indv=1:length(ider)
       vproc = outmode.(outmode.modes{indm}).var_proc{ider(indv)};
       switch vproc
       case 'gi';
          % sums over groups the deriving variable
          % for 'gi' option, assumes group dimension is the last
          derv = outmode.(outmode.modes{indm}).var_derv{ider(indv)};
          tempvar = outmode.(outmode.modes{indm}).(derv);
          idim = length(size(tempvar));
          % uses a nanmask to preserve nans - here adpts minimum nan number
          tmpnan = squeeze(all(isnan(tempvar),idim));
          tempvar1 = nansum(tempvar,idim);
          tempvar1(tmpnan) = nan;
          outn = outmode.(outmode.modes{indm}).var_outn{ider(indv)}; 
          outmode.(outmode.modes{indm}).(outn) = tempvar1;
       otherwise
          error(['processing mode ' vproc ' not valid']);
       end
    end
 end
      
%-----------------------------------------------------------------------------------------
% Adds in output structure
%-----------------------------------------------------------------------------------------
 boats.output = outmode;

%-----------------------------------------------------------------------------------------
% Save forcing arrays
%-----------------------------------------------------------------------------------------
 if strcmp(MAIN.sim_type,'h')
   boats.forcing_used.catchability = catchability_used;
   boats.forcing_used.price        = price_used;
   boats.forcing_used.cost_effort  = cost_effort_used;
   
   if strcmp(MAIN.reg_type,'omnT')
       boats.forcing_used.societenf    = FORC.societenf;
       boats.forcing_used.effEtarg     = FORC.effEtarg;
   end
 end

%-----------------------------------------------------------------------------------------
% Save time taken to run this script
%-----------------------------------------------------------------------------------------
  runtime_integrate       = toc;
  boats.runtime_integrate = runtime_integrate;
  
%-----------------------------------------------------------------------------------------
% Save regulation stats for each simulation
%-----------------------------------------------------------------------------------------
  if ~strcmp(MAIN.reg_type,'oa') && strcmp(MAIN.sim_type,'h')                  % Only if harvest is regulated (not 'oa')
    % Regulation stats (output)
    boats.output.stats.regulation_onset        = regulation_onset;
    boats.output.stats.indt_at_reg_onset       = indt_at_reg_onset;
    boats.output.stats.annual_harvest_max      = annual_harvest_max;
    
    % Regulation parameters
    boats.output.stats.reg_threshold      = ECON.reg_threshold;
    boats.output.stats.k_s                = ECON.k_s;
%     boats.output.stats.reg_onset          = ECON.reg_onset;
    
  end
 
 end
 
%**************************************************************************************************************
% END OF FUNCTION


%**************************************************************************************************************
% SUBFUNCTION
%**************************************************************************************************************
 function out = integrate_2d(tvar,area);
    tmp = area .* tvar;
    out = nansum(tmp(:));
 end
 
%**************************************************************************************************************
% END OF SCRIPT
