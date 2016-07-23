function definespatiotemporalvalues
%DEFINESPATIOTEMPORALVALUES Summary of this function goes here #TODO
%
%   Detailed explanation goes here #TODO
%
%   Last edit: Taryn Black, 24 June 2016

% ======================================================================= %

%   Global variables that are defined
    global RESOLUTION_XYZ GRIDPOINTS_X GRIDPOINTS_Y GRIDPOINTS_Z
    global TICKS_X TICKS_Y TICKS_Z TICKLABELS_X TICKLABELS_Y TICKLABELS_Z
    global nTIMESTEPS TIME_VECTOR
    
%   Global variables that are used
    global TIME_STOP TIME_START DT LENGTH HEIGHT WIDTH IMAX JMAX KMAX ...
        GHOSTCELLS_IJK NUMTICKLABELS_XYZ SCALE_FACTOR_XYZ
    
%%  Define number of timesteps to be processed and all individual timesteps
    nTIMESTEPS = (TIME_STOP-TIME_START)/DT;
    TIME_VECTOR = (TIME_START:nTIMESTEPS)*DT;
% ----------------------------------------------------------------------- %

    
%%  Define dimensional spatial resolution [meters per cell] and domain grid
    RESOLUTION_XYZ = [  LENGTH/(IMAX - GHOSTCELLS_IJK(1));
                        HEIGHT/(JMAX - GHOSTCELLS_IJK(2));
                        WIDTH/(KMAX - GHOSTCELLS_IJK(3))];
    GRIDPOINTS_X = RESOLUTION_XYZ(1):RESOLUTION_XYZ(1):LENGTH;
    GRIDPOINTS_Y = RESOLUTION_XYZ(2):RESOLUTION_XYZ(2):HEIGHT;
    GRIDPOINTS_Z = RESOLUTION_XYZ(3):RESOLUTION_XYZ(3):WIDTH;
% ----------------------------------------------------------------------- %


%%  Define tick locations and labels
    TICKS_X = linspace(0,LENGTH,NUMTICKLABELS_XYZ(1));
    TICKS_Y = linspace(0,HEIGHT,NUMTICKLABELS_XYZ(2));
    TICKS_Z = linspace(0,WIDTH,NUMTICKLABELS_XYZ(3));
    TICKLABELS_X = TICKS_X(2:end)/SCALE_FACTOR_XYZ(1);
    TICKLABELS_Y = TICKS_Y(2:end)/SCALE_FACTOR_XYZ(2);
    TICKLABELS_Z = TICKS_Z(2:end)/SCALE_FACTOR_XYZ(3);
% ----------------------------------------------------------------------- %

%% ===================================================================== %%
end

