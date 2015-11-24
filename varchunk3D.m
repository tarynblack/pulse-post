function [ varchunk ] = varchunk3D( fileID,importspec,IMAX,JMAX,KMAX,ghostcells )
%varchunk3D Summary of this function goes here
%   Detailed explanation goes here

    varchunk_t = textscan(fileID,importspec,IMAX*JMAX*KMAX,'Delimiter','\t');
    varchunk = cell2mat(varchunk_t);
    varchunk = reshape(varchunk,[JMAX IMAX KMAX]);
    varchunk = varchunk(ghostcells-1:end-(ghostcells/2),ghostcells-1:end-(ghostcells/2),ghostcells-1:end-(ghostcells/2));
    varchunk = permute(varchunk,[3 2 1]);

end

