function [ Stokes_S1,Stokes_S2,Stokes_S3 ] = calcStokes( RO_S1,RO_S2,RO_S3,...
    D_S1,D_S2,D_S3,vel_char,MU_G0,VENT_R,PULSE,FREQ )
%calcStokes calculates the Stokes number of each particle phase.
%   Detailed explanation goes here

    if strcmp(PULSE,'T') == 1
        Stokes_S1 = RO_S1*FREQ*(D_S1^2)/(18*MU_G0);
        Stokes_S2 = RO_S2*FREQ*(D_S2^2)/(18*MU_G0);
        Stokes_S3 = RO_S3*FREQ*(D_S3^2)/(18*MU_G0);
    elseif strcmp(PULSE,'F') == 1
        Stokes_S1 = RO_S1*(D_S1^2)*vel_char/(18*MU_G0*2*VENT_R);
        Stokes_S2 = RO_S2*(D_S2^2)*vel_char/(18*MU_G0*2*VENT_R);
        Stokes_S3 = RO_S3*(D_S3^2)*vel_char/(18*MU_G0*2*VENT_R);
    end


end

