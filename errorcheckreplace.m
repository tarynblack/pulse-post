%%% errorcheckreplace checks workspace variables for problems before
%%% processing, generates a warning or error if problems are found, and
%%% fixes the problem if possible. It also redefines user-defined values
%%% that need to be updated after loading all of the simulation parameters.
%
%   Last edit: Taryn Black, 24 June 2016

% ======================================================================= %

%%  Global variables that are defined
    global SLICE_AZIMUTH SLICE_ELEVATION

%%  Consistency between number of isosurfaces and their color and transparency
    [size_COLOR] = size(EPG_ISOSURFACE_COLOR);
    if size_COLOR(1) ~= length(EPG_ISOSURFACE_EDGE)
        error(['Error: Number of rows in EPG_ISOSURFACE_COLOR must ' ...
            'match the number of isosurfaces, defined by length of ' ...
            'EPG_ISOSURFACE_EDGE.'])
    elseif length(EPG_ISOSURFACE_TRANSPARENCY) ~= length(EPG_ISOSURFACE_EDGE)
        error(['Error: Length of EPG_ISOSURFACE_TRANSPARENCY must ' ...
            'match the number of isosurfaces, defined by length of ' ...
            'EPG_ISOSURFACE_EDGE.'])
    elseif size_COLOR(2) ~= 3
        error('Error: EPG_ISOSURFACE_COLOR must have 3 columns for RGB')
    end
    clear size_COLOR

    
%%  Replace maximum colorbar value for gas temperature with maximum initial 
%%  gas temperature in simulation.
    if length(COLORRANGE_GAS_TEMPERATURE) == 1 || ...
            isnan(COLORRANGE_GAS_TEMPERATURE(2))
        COLORRANGE_GAS_TEMPERATURE(2) = BC_TG;
    end
    
    
%%  Redefine mass flux altitudes as a vector of grid indices
    MASS_FLUX_ALTITUDES(MASS_FLUX_ALTITUDES == 0) = RESOLUTION_XYZ(2);
    MASS_FLUX_ALTITUDES(isnan(MASS_FLUX_ALTITUDES)) = ...
        round(JET_HEIGHT)*RESOLUTION_XYZ(2);
    MASS_FLUX_ALTITUDES = sort(MASS_FLUX_ALTITUDES);
    MASS_FLUX_LEGEND = strcat(strsplit(num2str(MASS_FLUX_ALTITUDES/SCALE_FACTOR_XYZ(2))),...
        {' '},cellstr(LABEL_UNIT_XYZ(2)));
    MASS_FLUX_ALTITUDES = MASS_FLUX_ALTITUDES./RESOLUTION_XYZ(2);
    
    
%%  Fix string syntax for RUN_ID to be used in figure labels
    RUN_NAME = strrep(RUN_ID,'_','\_');
    
    
%% Check format of SLICE values and determine appropriate viewing angles
    if SLICE_LOCATION_XYZ(1) == 0
        SLICE_LOCATION_XYZ(1) = NaN;
    end
    if SLICE_LOCATION_XYZ(2) == 0
        SLICE_LOCATION_XYZ(2) = NaN;
    end
    if SLICE_LOCATION_XYZ(3) == 0
        SLICE_LOCATION_XYZ(3) = NaN;
    end
    
    if isnan(SLICE_LOCATION_XYZ(1)) && isnan(SLICE_LOCATION_XYZ(2))
        SLICE_AZIMUTH = 0;
        SLICE_ELEVATION = 90;
    elseif isnan(SLICE_LOCATION_XYZ(2)) && isnan(SLICE_LOCATION_XYZ(3))
        SLICE_AZIMUTH = 90;
        SLICE_ELEVATION = 0;
    elseif isnan(SLICE_LOCATION_XYZ(1)) && isnan(SLICE_LOCATION_XYZ(3))
        SLICE_AZIMUTH = 0;
        SLICE_ELEVATION = 0;
    else [SLICE_AZIMUTH,SLICE_ELEVATION] = view(3);
    end