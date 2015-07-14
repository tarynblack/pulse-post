% testing load of binary data

clear all; close all; clc

cd rundata/378654

fid = fopen('EP_G_378654');
EP_G = fread(fid,'int16','a');

IMAX = 104;
JMAX = 104;
KMAX = 104;
timesteps = length(EP_G)/(IMAX*JMAX*KMAX);

t = 50;
EPG = reshape(EP_G((t-1)*IMAX*JMAX*KMAX+1:t*IMAX*JMAX*KMAX),[JMAX IMAX KMAX]);
EPG = EPG(3:end-2,3:end-2,3:end-2);
EPG = permute(EPG,[3 2 1]);

ylevel = 50;
EPGslice = squeeze(EPG(:,ylevel,:));
open EPGslice

cd ../..
