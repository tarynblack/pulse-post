%%% Converts multipage .tif file to set of image files.
%%% Taryn Black, 19 November 2015

%%% --------------------- USER-DEFINED PARAMETERS --------------------- %%%

% Define preferred image type extension (jpg, bmp, png...).
  newType = 'jpg';

% List .tif files to convert to new image type (include '.tif' extension).
  tifName = {'EPG_tsteps_mary_compare.tif';
             'Entr_tsteps_mary_compare.tif'};
         
% Define name prefix for output images, in same order as in tifName.
% Naming convention is: " varname_t###.newType "
  varname = {'EPG';
             'Entrainment'};
         
% Define path that new images will be saved under.
  savedir = '/Users/taryn/OneDrive/Documents/PULSE_SHARE/entrainment_comparison';
  
% Define time interval between tif frames (seconds).
  DT = 5;
  
%%% =================================================================== %%%

for i = 1:length(tifName)
    
    % Determine number of pages in .tif file
      tifInfo = imfinfo(tifName{i});
      numFrames = length(tifInfo);
      
    cd(savedir)
    
    % Step through tif frames and save as new image files.
      for frame = 1:numFrames
          tifFrame = imread(tifName{i},frame);
          imwrite(tifFrame,strcat(varname{i},'_t',num2str(frame*DT),newType));
      end
end