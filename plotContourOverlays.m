function plotContourOverlays( FIGURE,AXES,VAR,CONTOUR_VALUE,CONTOUR_COLOR )
%plotContourOverlays Summary of this function goes here
%   Detailed explanation goes here
%
%   Last edit: Taryn Black, 2 August 2016

%%  Global variables that are used
    global SLICE_LOCATION_XYZ IMAX JMAX KMAX GHOSTCELLS_IJK

%%  Replace NaN with ZERO in SLICE_LOCATION_XYZ
%   Otherwise contours plot as grids
    SLICE_LOCATION_XYZ(isnan(SLICE_LOCATION_XYZ)) = 0;
    
%%  Plot contour line overlay on existing axis/axes in figure.
    figure(FIGURE)
    
    switch length(AXES)
        
        case 1
            hold on
            
            hCONTOUR = contourslice(AXES,VAR,...
                SLICE_LOCATION_XYZ(:,1)*(IMAX - GHOSTCELLS_IJK(1)),...
                0*(KMAX - GHOSTCELLS_IJK(3)),...
                0*(JMAX - GHOSTCELLS_IJK(2)),...
                [CONTOUR_VALUE CONTOUR_VALUE]);
            set(hCONTOUR,'EdgeColor',CONTOUR_COLOR,'LineWidth',0.5);
            
        otherwise            
            for n = 1:length(AXES)
                hold on
                
                hCONTOUR = contourslice(AXES(n),VAR,...
                    SLICE_LOCATION_XYZ(:,1)*(IMAX - GHOSTCELLS_IJK(1)),...
                    SLICE_LOCATION_XYZ(:,3)*(KMAX - GHOSTCELLS_IJK(3)),...
                    SLICE_LOCATION_XYZ(:,2)*(JMAX - GHOSTCELLS_IJK(2)),...
                    [CONTOUR_VALUE CONTOUR_VALUE]);
                set(hCONTOUR,'EdgeColor',CONTOUR_COLOR,'LineWidth',0.5);
            end
    end


end

