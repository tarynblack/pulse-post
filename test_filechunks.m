clear all; close all; clc

IMAX = 104;
JMAX = 104;
KMAX = 104;
OUT_DT = 5;
timesteps = (600/OUT_DT)+1;

N = IMAX*JMAX*KMAX;
formatSpec = '%f';

fileID = fopen('EP_G_477394');

k = 0;
while ~feof(fileID)
    k = k+1;
    EPG_t = textscan(fileID,formatSpec,N,'Delimiter','\t');
    if k == 80
        EPG_80 = EPG_t;
    end
end