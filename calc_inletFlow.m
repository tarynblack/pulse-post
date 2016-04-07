function [ XG,vel_char,MFR,MASSFLUX,MFR_SOL,MASSFLUX_SOL,jetheight ] = ...
    calc_inletFlow( charEPG,MING,MAXG,PULSE,BC_EPG,BC_PG,BC_TG,Rgas,...
    RO_S1,RO_S2,RO_S3,NFR_S1,NFR_S2,NFR_S3,BC_TS1,BC_TS2,BC_TS3,VENT_R,g )
%calc_inletFlow calculates the gas mass fraction, choked velocity, and mass
%flow rate at the inlet (volcanic vent) for a simulation.
%   If flow is steady, the standard calculations of these values apply. 
%
%   If flow is unsteady, then in each calculation, the (steady) gas volume
%   fraction is replaced by a characteristic unsteady gas volume fraction
%   representing either the minimum, maximum, or average gas volume
%   fraction of a pulse. The user chooses which characteristic value to use
%   (charEPG). The resulting gass mass fraction, choked velocity, and mass
%   flow rate are then characteristic of the chosen gas volume fraction.
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
% Last edit: Taryn Black, 6 April 2016

    if strcmp(charEPG,'mingas') == 1
        EPGunst = MING;
    elseif strcmp(charEPG,'maxgas') == 1
        EPGunst = MAXG;
    elseif strcmp(charEPG,'avggas') == 1
        EPGunst = (2*MING*(1-MING)/pi)+MING;
    end

    if strcmp(PULSE,'T') == 1
        XG = (1 + ((1-EPGunst)*Rgas*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3)/(EPGunst*BC_PG)))^(-1);
        vel_char = sqrt((Rgas/XG)*(EPGunst*BC_TG +((1-EPGunst)*(BC_TS1*NFR_S1+BC_TS2*NFR_S2+BC_TS3*NFR_S3))))*(XG+((1-XG)*BC_PG/(Rgas*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3))));
        MFR = pi*(VENT_R^2)*vel_char*((EPGunst*BC_PG/(Rgas*BC_TG))+(1-EPGunst)*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3));
        MASSFLUX = MFR/(pi*(VENT_R^2));
        MFR_SOL = pi*(VENT_R^2)*vel_char*(1-EPGunst)*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3);
        MASSFLUX_SOL = MFR_SOL/(pi*(VENT_R^2));
    elseif strcmp(PULSE,'F') == 1
        XG = (1 + ((1-BC_EPG)*Rgas*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3)/(BC_EPG*BC_PG)))^(-1);
        vel_char = sqrt((Rgas/XG)*(BC_EPG*BC_TG +((1-BC_EPG)*(BC_TS1*NFR_S1+BC_TS2*NFR_S2+BC_TS3*NFR_S3))))*(XG+((1-XG)*BC_PG/(Rgas*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3))));
        MFR = pi*(VENT_R^2)*vel_char*((BC_EPG*BC_PG/(Rgas*BC_TG))+(1-BC_EPG)*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3));
        MASSFLUX = MFR/(pi*(VENT_R^2));
        MFR_SOL = pi*(VENT_R^2)*vel_char*(1-BC_EPG)*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3);
        MASSFLUX_SOL = MFR_SOL/(pi*(VENT_R^2));
    end
    
    jetheight = (vel_char^2)/(2*g);

end

