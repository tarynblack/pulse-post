function [ logEPS1,logEPS2,logEPS3 ] = ...
    calculateLogConcentrations( EPS1,EPS2,EPS3,shift_value )
%calculateLogConcentrations Summary of this function goes here
%   Detailed explanation goes here
%
%   Last edit: Taryn Black, 2 August 2016

logEPS1 = log10(EPS1 + shift_value);
logEPS2 = log10(EPS2 + shift_value);
logEPS3 = log10(EPS3 + shift_value);


end

