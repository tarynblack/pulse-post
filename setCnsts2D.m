function [IMAX,JMAX,HEIGHT,LENGTH,RO_S1,RO_S2,PULSE,FREQ] = setCnsts2D(run)
%setCnsts2D loads mfixconst_<run> and sets constant properties from an MFIX
%simulation specified by input argument <run>. The contents of mfixconst
%are specified in allocate_arrays.f in the run directory.
%   IMAX: number of cells in i (horiz) direction (plus ghost cells)
%   JMAX: number of cells in j (vert) direction (plus ghost cells)
%   HEIGHT: height of domain in meters
%   LENGTH: width of domain in meters
%   RO_S1: solid density of phase 1 particles
%   RO_S2: solid density of phase 2 particles
%   PULSE: string describes flow behavior. T = pulsing on; F = pulsing off.
%   FREQ: frequency of pulsing (valid only for PULSE = T). Note: FREQ will
%   return a value when onoff=F, but the value is unused in the
%   simulations.
%   Taryn Black, last edit 20 March 2015

cd(sprintf('%d',run))


% Read in second column of mfixconst_<run> as a cell array of strings
fid = fopen(sprintf('mfixconst_%d',run));
data = textscan(fid, '%*s %s');
fclose(fid);

% Convert to data structure (easier to save each cell as desired type)
Cnsts = cell2struct(data,{'data'},1);

% Save individual cells as either dbl or char
IMAX   = str2double(Cnsts.data(1))+2;
JMAX   = str2double(Cnsts.data(2))+2;
HEIGHT = str2double(Cnsts.data(3));
LENGTH = str2double(Cnsts.data(4));
RO_S1  = str2double(Cnsts.data(5));
RO_S2  = str2double(Cnsts.data(6));
PULSE  = char(Cnsts.data(7));
FREQ   = str2double(Cnsts.data(8));

cd ..

end

