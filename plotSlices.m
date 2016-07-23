function plotSlices( FIGURE,AXES,VAR,FACE_COLOR,...
    EDGE_COLOR,COLORBAR_RANGE,VARIABLE_NAME,t )
%plotSlices Summary of this function goes here #TODO
%   Detailed explanation goes here #TODO
%
%   Last edit: Taryn Black, 21 July 2016

%%  Global variables that are used
    global IMAX JMAX KMAX GHOSTCELLS_IJK SLICE_LOCATION_XYZ
    global PULSE TIME_VECTOR RUN_NAME FREQUENCY

%%  Plot single or multiple #TODO slices in single figure
    figure(FIGURE)
    
    switch length(AXES)
        
        case 1
            cla(AXES)
            hSLICE = slice(0.5:(IMAX - GHOSTCELLS_IJK(1) - 0.5),...
                           0.5:(KMAX - GHOSTCELLS_IJK(3) - 0.5),...
                           0.5:(JMAX - GHOSTCELLS_IJK(2) - 0.5),...
                           VAR,...
                           SLICE_LOCATION_XYZ(1)*(IMAX - GHOSTCELLS_IJK(1)),...
                           SLICE_LOCATION_XYZ(3)*(KMAX - GHOSTCELLS_IJK(3)),...
                           SLICE_LOCATION_XYZ(2)*(JMAX - GHOSTCELLS_IJK(2)));
            set(hSLICE,'FaceColor',FACE_COLOR,...
                       'EdgeColor',EDGE_COLOR);
    
            caxis(AXES,COLORBAR_RANGE)
    
            FIGURE_TITLE = pulsetitle(VARIABLE_NAME,PULSE,TIME_VECTOR,t,RUN_NAME,FREQUENCY);
            title(FIGURE_TITLE,'FontWeight','bold');
            
        otherwise
            POSITION = zeros(length(AXES),4);
            for n = 1:length(AXES)
                cla(AXES(n))
                hold on
                hSLICE = slice(AXES(n),0.5:(IMAX - GHOSTCELLS_IJK(1) - 0.5),...
                                    0.5:(KMAX - GHOSTCELLS_IJK(3) - 0.5),...
                                    0.5:(JMAX - GHOSTCELLS_IJK(2) - 0.5),...
                                    VAR{n},...
                                    SLICE_LOCATION_XYZ(1)*(IMAX - GHOSTCELLS_IJK(1)),...
                                    SLICE_LOCATION_XYZ(3)*(KMAX - GHOSTCELLS_IJK(3)),...
                                    SLICE_LOCATION_XYZ(2)*(JMAX - GHOSTCELLS_IJK(2)));
                set(hSLICE,'FaceColor',FACE_COLOR,...
                           'EdgeColor',EDGE_COLOR);
                POSITION(n,:) = get(AXES(n),'position');
                POSITION(n,3:4) = POSITION(1,3:4);
                set(AXES(n),'position',POSITION(n,:));
                caxis(AXES(n),COLORBAR_RANGE)
                FIGURE_TITLE = pulsetitle(VARIABLE_NAME{n},PULSE,TIME_VECTOR,t,RUN_NAME,FREQUENCY);
                title(AXES(n),FIGURE_TITLE,'FontWeight','bold');
            end
            
    end
                  

end

