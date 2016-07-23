function [ T_G ] = calculateAdjustedTemperature( T_G )
%calculateAdjustedTemperature Summary of this function goes here #TODO
%   Detailed explanation goes here #TODO
%   Apply atmospheric correction (equation of state for gas)...
%
%   Last edit: Taryn Black, 21 July 2016

%%  Global variables that are used
    global ATMOS TROPO JMAX GHOSTCELLS_IJK GRIDPOINTS_Y
    
%%  Atmospheric correction for temperature
    if strcmp(ATMOS,'T') == 1
        for j = 1:(JMAX - GHOSTCELLS_IJK)
            if GRIDPOINTS_Y(j) <= TROPO
                T_G(:,:,j) = T_G(:,:,j) - 0.0098*GRIDPOINTS_Y(j);
            elseif GRIDPOINTS_Y(j) > TROPO
                T_G(:,:,j) = T_G(:,:,j) - 0.0098*TROPO + 0.001*GRIDPOINTS_Y(j);
            end
        end
    end


end

