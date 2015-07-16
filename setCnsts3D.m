function [ IMAX,JMAX,KMAX,LENGTH,HEIGHT,WIDTH,RO_S1,RO_S2,PULSE,FREQ,...
    DT ] = setCnsts3D( run,dir )
%setCnsts3D loads the file mfixconst_<run> and sets constant parameters for
%post-processing from an MFiX simulation specified by <run>. The contents
%of mfixconst are specified in allocate_arrays.f in the simulation run
%directory.
%   IMAX:   number of cells in x direction
%   JMAX:   number of cells in y direction
%   KMAX:   number of cells in z direction
%       *** NOTE: this function adds 4 domain ghost cells to each *MAX
%   LENGTH: size of domain in x direction (meters)
%   HEIGHT: size of domain in y direction (meters)
%   WIDTH:  size of domain in z direction (meters)
%   RO_S#:  solid density of phase# particles (kg/m3)
%   PULSE:  describes whether flow is pulsing (T) or steady (F)
%   FREQ:   frequency of pulsing, for PULSE=T
%       *** NOTE: FREQ returns a value for PULSE=F that is not used.
%   DT:     time interval between each data write in the simulation
    
    cd(dir)
    
    % Read in second column of mfixconst (containing values) as a cell
    % array of strings
    fid = fopen(sprintf('mfixconst_%d',run));
    data = textscan(fid, '%*s %s');
    fclose(fid);
    
    % Convert to data structure to save each cell as desired type
    Cnsts = cell2struct(data,{'data'},1);
    ghostcells = 4;
    
    IMAX   = str2double(Cnsts.data(1)) + ghostcells;
    JMAX   = str2double(Cnsts.data(2)) + ghostcells;
    KMAX   = str2double(Cnsts.data(3)) + ghostcells;
    LENGTH = str2double(Cnsts.data(4));
    HEIGHT = str2double(Cnsts.data(5));
    WIDTH  = str2double(Cnsts.data(6));
    RO_S1  = str2double(Cnsts.data(7));
    RO_S2  = str2double(Cnsts.data(8));
    PULSE  = char(Cnsts.data(9));
    FREQ   = str2double(Cnsts.data(10));
    DT     = str2double(Cnsts.data(11));

    
end