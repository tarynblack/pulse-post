function [ Stokes_S1,Stokes_S2,Stokes_S3 ] = ...
    computeparticlestokes( CHARACTERISTIC_VELOCITY)
%COMPUTEPARTICLESTOKES calculates the Stokes number of each particle phase.
%   Detailed explanation goes here #TODO

%   Global variables that are used
    global RO_S1 RO_S2 RO_S3 D_S1 D_S2 D_S3 MU_G0 VENT_R PULSE FREQUENCY SAVEPATH 
    
if strcmp(PULSE,'T') == 1
    Stokes_S1 = RO_S1*FREQUENCY*(D_S1^2)/(18*MU_G0);
    Stokes_S2 = RO_S2*FREQUENCY*(D_S2^2)/(18*MU_G0);
    Stokes_S3 = RO_S3*FREQUENCY*(D_S3^2)/(18*MU_G0);
elseif strcmp(PULSE,'F') == 1
    Stokes_S1 = RO_S1*(D_S1^2)*CHARACTERISTIC_VELOCITY/(18*MU_G0*2*VENT_R);
    Stokes_S2 = RO_S2*(D_S2^2)*CHARACTERISTIC_VELOCITY/(18*MU_G0*2*VENT_R);
    Stokes_S3 = RO_S3*(D_S3^2)*CHARACTERISTIC_VELOCITY/(18*MU_G0*2*VENT_R);
end
    
dlmwrite(fullfile(SAVEPATH,'ParticleStokesNumbers.txt'),...
    [Stokes_S1 Stokes_S2 Stokes_S3],'delimiter','\t');

end
