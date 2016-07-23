function saveImageFile( FIGURE,VARIABLE_NAME,t )
%saveImageFile Summary of this function goes here #TODO
%   Detailed explanation goes here #TODO
%
%   Last edit: Taryn Black, 21 July 2016

%%  Global variables that are used
    global IMAGE_SAVE_TYPE SAVEPATH TIME_VECTOR RUN_NAME

%%  Save figure as frame of multipage TIF file if user specifies TIF 
%%  filetype, or save figure as single image if user specifies other type.
    if strcmp(IMAGE_SAVE_TYPE,'tif') == 1 || strcmp(IMAGE_SAVE_TYPE,'tiff') == 1
        IMAGE = getframe(FIGURE);
        imwrite(frame2im(IMAGE),fullfile(SAVEPATH,...
            sprintf('%s_allTimeSteps_%s.tif',VARIABLE_NAME,RUN_NAME)),...
            'tif','WriteMode','append');
    else
        saveas(FIGURE,fullfile(SAVEPATH,sprintf('%s_time%03d_%s.%s',...
            VARIABLE_NAME,TIME_VECTOR(t),RUN_NAME,IMAGE_SAVE_TYPE)));
    end
    

end

