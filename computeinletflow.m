function [ AVERAGE_EPG,GAS_MASS_FRACTION,CHARACTERISTIC_VELOCITY,...
    MASS_FLOW_RATE,MASS_FLUX,SOLID_MASS_FLOW_RATE,SOLID_MASS_FLUX,...
    JET_HEIGHT] = computeinletflow
%COMPUTEINLETFLOW calculates the gas mass fraction, choked velocity, and
%mass flow rate at the inlet (volcanic vent) for a simulation.
%
%   If flow is steady, the standard calculations of these values apply. 
%
%   If flow is unsteady, then in each calculation, the (steady) gas volume
%   fraction is replaced by a characteristic unsteady gas volume fraction
%   representing either the minimum, maximum, or average gas volume
%   fraction of a pulse. The user chooses which characteristic value to use
%   (CHARACTERISTIC_EPG). The resulting gass mass fraction, choked
%   velocity, and mass flow rate are then characteristic of the chosen gas
%   volume fraction.
%   
%   If average gas volume fraction is selected, this function uses the
%   result of the integral of the pulse over the pulse period, where the
%   pulse has the form: 
%     EPG(t) = MING*(1-MING)*abs(sin(2*pi*FREQ*t)) + MING 
%   to calculate the average gas volume fraction. Note that this pulse form
%   is fully sinusoidal and not clipped at values above MAXG as in the
%   actual simulation. For large MAXG (close to 1) the difference should
%   not be significant.
%
% Last edit: Taryn Black, 24 June 2016

% add comments describing calculations #TODO

%   Global variables that are used
    global CHARACTERISTIC_EPG MINIMUM_EPG MAXIMUM_EPG PULSE BC_EPG ...
        BC_PG BC_TG R_GAS RO_S1 RO_S2 RO_S3 NFR_S1 NFR_S2 NFR_S3 BC_TS1 ...
        BC_TS2 BC_TS3 VENT_R g RESOLUTION_XYZ

    BETA = asin((MAXIMUM_EPG-MINIMUM_EPG)/(MINIMUM_EPG*(1-MINIMUM_EPG)));
    AVERAGE_EPG = 1 + ((2/pi)*(1-MINIMUM_EPG)*(MINIMUM_EPG - (MINIMUM_EPG*cos(BETA)) - BETA));

    if strcmp(CHARACTERISTIC_EPG,'minimumEPG') == 1
        UNSTEADY_EPG = MINIMUM_EPG;
    elseif strcmp(CHARACTERISTIC_EPG,'maximumEPG') == 1
        UNSTEADY_EPG = MAXIMUM_EPG;
    elseif strcmp(CHARACTERISTIC_EPG,'averageEPG') == 1
        UNSTEADY_EPG = AVERAGE_EPG;
    end

    if strcmp(PULSE,'T') == 1
        GAS_MASS_FRACTION = (1 + ((1-UNSTEADY_EPG)*R_GAS*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3)/(UNSTEADY_EPG*BC_PG)))^(-1);
        CHARACTERISTIC_VELOCITY = sqrt((R_GAS/GAS_MASS_FRACTION)*(UNSTEADY_EPG*BC_TG +((1-UNSTEADY_EPG)*(BC_TS1*NFR_S1+BC_TS2*NFR_S2+BC_TS3*NFR_S3))))*(GAS_MASS_FRACTION+((1-GAS_MASS_FRACTION)*BC_PG/(R_GAS*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3))));
        MASS_FLOW_RATE = pi*(VENT_R^2)*CHARACTERISTIC_VELOCITY*((UNSTEADY_EPG*BC_PG/(R_GAS*BC_TG))+(1-UNSTEADY_EPG)*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3));
        MASS_FLUX = MASS_FLOW_RATE/(pi*(VENT_R^2));
        SOLID_MASS_FLOW_RATE = pi*(VENT_R^2)*CHARACTERISTIC_VELOCITY*(1-UNSTEADY_EPG)*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3);
        SOLID_MASS_FLUX = SOLID_MASS_FLOW_RATE/(pi*(VENT_R^2));
    elseif strcmp(PULSE,'F') == 1
        GAS_MASS_FRACTION = (1 + ((1-BC_EPG)*R_GAS*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3)/(BC_EPG*BC_PG)))^(-1);
        CHARACTERISTIC_VELOCITY = sqrt((R_GAS/GAS_MASS_FRACTION)*(BC_EPG*BC_TG +((1-BC_EPG)*(BC_TS1*NFR_S1+BC_TS2*NFR_S2+BC_TS3*NFR_S3))))*(GAS_MASS_FRACTION+((1-GAS_MASS_FRACTION)*BC_PG/(R_GAS*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3))));
        MASS_FLOW_RATE = pi*(VENT_R^2)*CHARACTERISTIC_VELOCITY*((BC_EPG*BC_PG/(R_GAS*BC_TG))+(1-BC_EPG)*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3));
        MASS_FLUX = MASS_FLOW_RATE/(pi*(VENT_R^2));
        SOLID_MASS_FLOW_RATE = pi*(VENT_R^2)*CHARACTERISTIC_VELOCITY*(1-BC_EPG)*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3);
        SOLID_MASS_FLUX = SOLID_MASS_FLOW_RATE/(pi*(VENT_R^2));
    end
    
    JET_HEIGHT = (CHARACTERISTIC_VELOCITY^2)/(2*g)/RESOLUTION_XYZ(2);

end

