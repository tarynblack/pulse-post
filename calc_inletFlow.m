function [ XG,vel_char,MFR ] = calc_inletFlow( charEPG,MING,MAXG,PULSE,...
    BC_EPG_st,BC_PG,BC_TG,Rgas,RO_S1,RO_S2,RO_S3,NFR_S1,NFR_S2,NFR_S3,...
    BC_TS1,BC_TS2,BC_TS3,VENT_R )
%calc_inletFlow calculates the gas mass fraction, choked velocity, and mass
%flow rate at the inlet (volcanic vent) for a simulation.
%   Detailed explanation goes here
% Last edit: Taryn Black, 14 November 2015

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
    elseif strcmp(PULSE,'F') == 1
        XG = (1 + ((1-BC_EPG_st)*Rgas*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3)/(BC_EPG_st*BC_PG)))^(-1);
        vel_char = sqrt((Rgas/XG)*(BC_EPG_st*BC_TG +((1-BC_EPG_st)*(BC_TS1*NFR_S1+BC_TS2*NFR_S2+BC_TS3*NFR_S3))))*(XG+((1-XG)*BC_PG/(Rgas*BC_TG*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3))));
        MFR = pi*(VENT_R^2)*vel_char*((BC_EPG_st*BC_PG/(Rgas*BC_TG))+(1-BC_EPG_st)*(RO_S1*NFR_S1+RO_S2*NFR_S2+RO_S3*NFR_S3));
    end

end

