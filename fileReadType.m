function [ fID ] = fileReadType( fname,readtype,t,runpath,postpath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cd(runpath)

if readtype == 1
    fID = fopen(sprintf('%s_t%02d.txt',fname{readtype},t));
elseif readtype == 2
    fID = fopen(sprintf('%s',fname{readtype}));
end

cd(postpath)

end

