function [ SOLID_MASS_FLUX_V,...
           NET_SOLID_MASS_FLUX_V,...
           NET_SOLID_MASS_FLUX_V_USER_ALTITUDES,...
           LOG_SOLID_MASS_FLUX_V ] ...
           = calculateMassFlux( EPS1,EPS2,EPS3,V_S1,V_S2,V_S3 )
%calculateMassFlux Summary of this function goes here #TODO
%   Detailed explanation goes here #TODO
%
%   Last edit: Taryn Black, 22 July 2016

%%  Global variables that are used
    global RO_S1 RO_S2 RO_S3
    global MASS_FLUX_ALTITUDES

%%  Calculate vertical solid mass flux (total and net)
    SOLID_MASS_FLUX_V = RO_S1*EPS1.*V_S1 + RO_S2*EPS2.*V_S2 + RO_S3*EPS3.*V_S3;
    NET_SOLID_MASS_FLUX_V = squeeze(sum(sum(SOLID_MASS_FLUX_V)));
    
%%  Select net mass flux values at altitudes specified by user
    NET_SOLID_MASS_FLUX_V_USER_ALTITUDES = ...
        NET_SOLID_MASS_FLUX_V(MASS_FLUX_ALTITUDES);
      
%%  Separate positive and negative mass fluxes for logarithmic plotting
    LOG_SOLID_MASS_FLUX_V = SOLID_MASS_FLUX_V;
    LOG_SOLID_MASS_FLUX_V(LOG_SOLID_MASS_FLUX_V > 0) = ...
        log10(LOG_SOLID_MASS_FLUX_V(LOG_SOLID_MASS_FLUX_V > 0));
    LOG_SOLID_MASS_FLUX_V(LOG_SOLID_MASS_FLUX_V < 0) = ...
        -log10(abs(LOG_SOLID_MASS_FLUX_V(LOG_SOLID_MASS_FLUX_V < 0)));


end

