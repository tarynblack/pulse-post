function [ timeslice ] = timeslice3D( full_var,t,IMAX,JMAX,KMAX,ghostcells )
%timeslice3D returns a 3D matrix of values for an MFiX variable at a single
%timestep.
%   timeslice3D takes an entire MFiX space-time vector <full_var>, reshapes
%   it into a 3D matrix <timeslice> for a single timestep <t>, removes the
%   <ghostcells> and permutes it to the correct Cartesian coordinates.

    timeslice = reshape(full_var((t-1)*IMAX*JMAX*KMAX+1:t*IMAX*JMAX*KMAX),[JMAX IMAX KMAX]);
    timeslice = timeslice(ghostcells-1:end-(ghostcells/2),ghostcells-1:end-(ghostcells/2),ghostcells-1:end-(ghostcells/2));
    timeslice = permute(timeslice,[3 2 1]);

end

