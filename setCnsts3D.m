function [ IMAX,JMAX,KMAX,LENGTH,HEIGHT,WIDTH,RO_S1,RO_S2,RO_S3,NFR_S1,...
    NFR_S2,NFR_S3,PULSE,FREQ,MING,MAXG,VENT_R,DT,TSTOP,ATMOS,TROPO,...
    BC_EPG,BC_PG,BC_TG,BC_TS1,BC_TS2,BC_TS3 ] = setCnsts3D( run,dir,ghostcells )
% setCnsts3D loads the file mfixconst_<run> and sets constant parameters for
% post-processing from an MFiX simulation specified by <run>. The contents
% of mfixconst are specified in allocate_arrays.f in the simulation run
% directory.
%   IMAX:   number of cells in x direction
%   JMAX:   number of cells in y direction
%   KMAX:   number of cells in z direction
%       *** NOTE: this function adds 4 domain ghost cells to each *MAX
%   LENGTH: size of domain in x direction (meters)
%   HEIGHT: size of domain in y direction (meters)
%   WIDTH:  size of domain in z direction (meters)
%   RO_S#:  solid density of phase# particles (kg/m3)
%   NFR_S#: number fraction of phase# particles (sum of NFR_S# == 1)
%   PULSE:  describes whether flow is pulsing (T) or steady (F)
%   FREQ:   frequency of pulsing, for PULSE=T
%       *** NOTE: FREQ returns a value for PULSE=F that is not used.
%   MING:   minimum gas volume fraction (unsteady flow)
%   MAXG:   maximum gas volume fraction (unsteady flow)
%   VENT_R: radius of vent [m]
%   DT:     time interval between each data write in the simulation
%   TSTOP:  end time of simulation [s]
%   ATMOS:  describes whether simulation runs with atmospheric temperature,
%           density, etc. conditions (T) or not (F)
%   TROPO:  altitude of tropopause [m]
%   BC_EPG: inlet gas volume fraction (constant for steady case)
%   BC_PG:  inlet gas pressure [Pa, or N/m2]
%   BC_TG:  inlet gas temperature [K]
%   BC_TS#: inlet temperature of phase# particles [K]
%
% Last edit: Taryn Black, 17 November 2015
    
    cd(dir)
    
    % Read in second column of mfixconst (containing values) as a cell
    % array of strings
    fid = fopen(sprintf('mfixconst_%d',run));
    data = textscan(fid, '%*s %s');
    fclose(fid);
    
    % Convert to data structure to save each cell as desired type
    Cnsts = cell2struct(data,{'data'},1);
    
    IMAX   = str2double(Cnsts.data(1)) + ghostcells;
    JMAX   = str2double(Cnsts.data(2)) + ghostcells;
    KMAX   = str2double(Cnsts.data(3)) + ghostcells;
    LENGTH = str2double(Cnsts.data(4));
    HEIGHT = str2double(Cnsts.data(5));
    WIDTH  = str2double(Cnsts.data(6));
    RO_S1  = str2double(Cnsts.data(7));
    RO_S2  = str2double(Cnsts.data(8));
    RO_S3  = str2double(Cnsts.data(9));
    NFR_S1 = str2double(Cnsts.data(10));
    NFR_S2 = str2double(Cnsts.data(11));
    NFR_S3 = str2double(Cnsts.data(12));
    PULSE  = char(Cnsts.data(13));
    FREQ   = str2double(Cnsts.data(14));
    MING   = str2double(Cnsts.data(15));
    MAXG   = str2double(Cnsts.data(16));
    VENT_R = str2double(Cnsts.data(17));
    DT     = str2double(Cnsts.data(18));
    TSTOP  = str2double(Cnsts.data(19));
    ATMOS  = char(Cnsts.data(20));
    TROPO  = str2double(Cnsts.data(21));
    BC_EPG = str2double(Cnsts.data(22));
    BC_PG  = str2double(Cnsts.data(23));
    BC_TG  = str2double(Cnsts.data(24));
    BC_TS1 = str2double(Cnsts.data(25));
    BC_TS2 = str2double(Cnsts.data(26));
    BC_TS3 = str2double(Cnsts.data(27));

    
end