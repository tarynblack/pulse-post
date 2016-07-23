function plotIsosurfaces( FIGURE,AXES,VAR,...
    ISOSURFACE_VALUE,ISOSURFACE_COLOR,ISOSURFACE_EDGE_COLOR,...
    ISOSURFACE_TRANSPARENCY,ISOSURFACE_LIGHTING,VARIABLE_NAME,t )
%plotIsosurfaces Summary of this function goes here #TODO
%   Detailed explanation goes here #TODO
%
%   Last edit: Taryn Black, 21 July 2016

%%  Global variables that are used
    global PULSE TIME_VECTOR RUN_NAME FREQUENCY

    
%%  Plot single or multiple isosurfaces in single axes
    figure(FIGURE)
    cla(AXES)
    
    for n = 1:length(ISOSURFACE_VALUE)
        surf(n) = patch(isosurface(VAR,ISOSURFACE_VALUE(n)));        
        set(surf(n),...
            'FaceColor',ISOSURFACE_COLOR(n,:),...
            'EdgeColor',ISOSURFACE_EDGE_COLOR,...
            'FaceAlpha',ISOSURFACE_TRANSPARENCY(n));
        hold on
    end
    
%%  Add figure lighting if LIGHTING is 'on'
    if strcmp(ISOSURFACE_LIGHTING,'on') == 1
        camlight('right');
        camlight('left');
        lighting gouraud;
    end

%%  Add figure title
    FIGURE_TITLE = pulsetitle(VARIABLE_NAME,PULSE,TIME_VECTOR,t,RUN_NAME,FREQUENCY);
    title(FIGURE_TITLE,'FontWeight','bold');
    set(FIGURE,'Visible','on');

end

