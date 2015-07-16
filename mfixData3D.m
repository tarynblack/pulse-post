function [ EP_G,...
            P_G,...
            T_G,T_S1,T_S2,...
            U_G,U_S1,U_S2,...
            V_G,V_S1,V_S2,...
            W_G,W_S1,W_S2,...
            RO_G,ROP_S1,ROP_S2,...
            X_G2 ] = mfixData3D( run,dir,onE,onP,onT,onV,onR,onX )
% mfixData3D loads output data from an MFiX simulation specified by <run>.
%   If on* == 1, load the data vector. If on* == 0, create an empty vector.
%
%   Vectors represent:
%   EP_G:   gas volume fraction (onE)
%   P_G:    gas pressure (onP)
%   T_*:    temperature (*: G - gas; S1 - solid1; S2 - solid2) (onT)
%   U_*:    velocity in the x direction (onV)
%   V_*:    velocity in the y direction (onV)
%   W_*:    velocity in the z direction (onV)
%   RO_G:   effective gas density (onR)
%   ROP_S#: effective phase# density (onR)
%   X_G2:   water vapor mass ratio (onX)
%
%   Last edit: Taryn Black, 16 July 2015

    cd(dir)
    
    if onE == 1
        EP_G = importdata(sprintf('EP_G_%d',run));
    elseif onE == 0
        EP_G = [];
    end
    
    if onP == 1
        P_G = importdata(sprintf('P_G_%d',run));
    elseif onP == 0
        P_G = [];
    end
    
    if onT == 1
        T_G  = importdata(sprintf('T_G_%d',run));
        T_S1 = importdata(sprintf('T_S1_%d',run));
        T_S2 = importdata(sprintf('T_S2_%d',run));
    elseif onT == 0
        T_G  = [];
        T_S1 = [];
        T_S2 = [];
    end
          
    if onV == 1
        U_G  = importdata(sprintf('U_G_%d',run));
        U_S1 = importdata(sprintf('U_S1_%d',run));
        U_S2 = importdata(sprintf('U_S2_%d',run));
        V_G  = importdata(sprintf('V_G_%d',run));
        V_S1 = importdata(sprintf('V_S1_%d',run));
        V_S2 = importdata(sprintf('V_S2_%d',run));
        W_G  = importdata(sprintf('W_G_%d',run));
        W_S1 = importdata(sprintf('W_S1_%d',run));
        W_S2 = importdata(sprintf('W_S2_%d',run));
    elseif onV == 0
        U_G  = [];
        U_S1 = [];
        U_S2 = [];
        V_G  = [];
        V_S1 = [];
        V_S2 = [];
        W_G  = [];
        W_S1 = [];
        W_S2 = [];
    end
        
    if onR == 1
        RO_G   = importdata(sprintf('RO_G_%d',run));
        ROP_S1 = importdata(sprintf('ROP_S1_%d',run));
        ROP_S2 = importdata(sprintf('ROP_S2_%d',run));
    elseif onR == 0
        RO_G   = [];
        ROP_S1 = [];
        ROP_S2 = [];
    end
        
    if onX == 1
        X_G2  = importdata(sprintf('X_G2_%d',run));
    elseif onX == 0
        X_G2 = [];
    end


end

