function [EP_G,RO_G,ROP_S1,ROP_S2,V_G,U_G,V_S1,V_S2,T_G,P_G] = mfixData(run)
% mfixData2D loads output data from an MFIX simulation specified by <run>
%   EP_G: gas volume fraction
%   RO_G: effective gas density
%   ROP_S1: effective solid(1) density
%   ROP_S2: effective solid(2) density
%   V_G: vertical gas velocity
%   V_S1: vertical solid(1) velocity
%   V_S2: vertical solid(2) velocity
%   T_G: gas temperature
%   P_G: gas pressure
%   Taryn Black, last edit 20 March 2015

cd(sprintf('%d',run))

    EP_G = importdata(sprintf('EP_G_%d',run));
    RO_G = importdata(sprintf('RO_G_%d',run));
    ROP_S1 = importdata(sprintf('ROP_S1_%d',run));
    ROP_S2 = importdata(sprintf('ROP_S2_%d',run));
    V_G = importdata(sprintf('V_G_%d',run));
    U_G = importdata(sprintf('U_G_%d',run));
    V_S1 = importdata(sprintf('V_S1_%d',run));
    V_S2 = importdata(sprintf('V_S2_%d',run));
    T_G = importdata(sprintf('T_G_%d',run));
    P_G = importdata(sprintf('P_G_%d',run));
    
cd ..

end

