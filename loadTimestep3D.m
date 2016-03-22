function [ varTimestep ] = loadTimestep3D( fileID,importspec,readtype,...
    IMAX,JMAX,KMAX,ghostcells )
%loadTimestep3D Summary of this function goes here
%   readtype == 1 :  file was originally written in binary, and was
%       converted to ascii using convert_Pulse.f90
%   readtype == 2 :  file was originally written in ascii and exists as a
%       single 'super-file' instead of pre-separated timesteps.

    if readtype == 1
        varTimestep_t = textscan(fileID,importspec,IMAX*JMAX*KMAX,'Delimiter','\t');
    elseif readtype == 2
        varTimestep_t = textscan(fileID,'%f',IMAX*JMAX*KMAX,'Delimiter','\t');
    end
    varTimestep = cell2mat(varTimestep_t);
    varTimestep = reshape(varTimestep,[JMAX IMAX KMAX]);
    varTimestep = varTimestep(ghostcells-1:end-(ghostcells/2),ghostcells-1:end-(ghostcells/2),ghostcells-1:end-(ghostcells/2));
    varTimestep = permute(varTimestep,[3 2 1]);
   
end

