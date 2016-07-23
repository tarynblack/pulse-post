function [ FLOW_VELOCITY,FLOW_VORTICITY_X,FLOW_VORTICITY_Y ] ...
    = calculateFlowVelocityVorticity( U_G,V_G,W_G )
%calculateFlowVelocityVorticity Summary of this function goes here #TODO
%   Detailed explanation goes here #TODO
%
%   Last edit: Taryn Black, 21 June 2016

FLOW_VELOCITY = sqrt(U_G.^2 + V_G.^2 + W_G.^2);

[FLOW_VORTICITY_X,FLOW_VORTICITY_Y] = curl(U_G,W_G,V_G);

end

