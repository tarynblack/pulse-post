function [ FIGURE,AXES,COLORBAR ] = createHandles( AZIMUTH,ELEVATION,NAME,...
    OUTER_POSITION,COLORMAP,CAXIS,nAXES )
%createHandles Summary of this function goes here #TODO
%   Detailed explanation goes here #TODO
%
%   Last edit: Taryn Black, 21 July 2016

%%  Argument management
    if nargin < 7 || nAXES < 1
        nAXES = 1;
    end

%%  Global variables that are used
    global IMAX JMAX KMAX GHOSTCELLS_IJK TICKS_X TICKS_Y TICKS_Z ...
        RESOLUTION_XYZ TICKLABELS_X TICKLABELS_Y TICKLABELS_Z LABEL_UNIT_XYZ

%%  Create figure, axes, and colorbar handles and set their properties
    switch nAXES
        
    %   Only one axes (single graph) in figure
        case 1
            FIGURE = figure('Name',NAME,'OuterPosition',OUTER_POSITION);
            AXES = axes('Parent',FIGURE);
            hold on
            axis(AXES,'equal',[0,IMAX - GHOSTCELLS_IJK(1),...
                               0,KMAX - GHOSTCELLS_IJK(3),...
                               0,JMAX - GHOSTCELLS_IJK(2)]);
            set(AXES,'XTick',TICKS_X(2:end)/RESOLUTION_XYZ(1),...
                        'XTickLabel',TICKLABELS_X,...
                     'YTick',TICKS_Z(2:end)/RESOLUTION_XYZ(3),...
                        'YTickLabel',TICKLABELS_Z,...
                     'ZTick',TICKS_Y(2:end)/RESOLUTION_XYZ(2),...
                        'ZTickLabel',TICKLABELS_Y); 
            grid(AXES,'on');AXES.Layer = 'top';
            view(AXES,AZIMUTH,ELEVATION);
            xlabel(AXES,sprintf('\\bf Distance_x (%s)',LABEL_UNIT_XYZ{1}));
            ylabel(AXES,sprintf('\\bf Distance_z (%s)',LABEL_UNIT_XYZ{3}));
            zlabel(AXES,sprintf('\\bf Altitude (%s)',LABEL_UNIT_XYZ{2}));
            
            colormap(AXES,COLORMAP);
            COLORBAR = colorbar(AXES(nAXES),'AxisLocation','in','FontSize',12);
            caxis(CAXIS);
        
            
    %   Multiple axes in figure (utilizes subtightplot)
        otherwise
            POSITION = zeros(nAXES,4);
            
            FIGURE = figure('Name',NAME,'OuterPosition',OUTER_POSITION);
            AXES = repmat(axes('Parent',FIGURE),1,nAXES);
            hold on
            axis(AXES,'equal',[0,IMAX - GHOSTCELLS_IJK(1),...
                               0,KMAX - GHOSTCELLS_IJK(3),...
                               0,JMAX - GHOSTCELLS_IJK(2)]);
            set(AXES,'XTick',TICKS_X(2:end)/RESOLUTION_XYZ(1),...
                        'XTickLabel',TICKLABELS_X,...
                     'YTick',TICKS_Z(2:end)/RESOLUTION_XYZ(3),...
                        'YTickLabel',TICKLABELS_Z,...
                     'ZTick',TICKS_Y(2:end)/RESOLUTION_XYZ(2),...
                        'ZTickLabel',TICKLABELS_Y);
            for n = 1:nAXES
                hold on
                AXES(n) = subtightplot(1,nAXES,n,[0.01 0.01],0.15,0.15);
                grid(AXES(n),'on');AXES(n).Layer = 'top';
                view(AXES(n),AZIMUTH,ELEVATION);
                xlabel(AXES(n),sprintf('\\bf Distance_x (%s)',LABEL_UNIT_XYZ{1}));
                ylabel(AXES(n),sprintf('\\bf Distance_z (%s)',LABEL_UNIT_XYZ{3}));
                zlabel(AXES(n),sprintf('\\bf Altitude (%s)',LABEL_UNIT_XYZ{2}));                
                colormap(AXES(n),COLORMAP);
                POSITION(n,:) = get(AXES(n),'position');
                POSITION(n,3:4) = POSITION(1,3:4);
                set(AXES(n),'position',POSITION(n,:));
            end
            COLORBAR = colorbar(AXES(nAXES),'AxisLocation','in','FontSize',12);
            caxis(CAXIS);
            
    end
     

end

