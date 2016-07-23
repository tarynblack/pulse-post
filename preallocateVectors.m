function [ surf_EPG,...
    AVG_EEE_COEFF,STD_EEE_COEFF,...
    AVG_ENTRAINMENT,STD_ENTRAINMENT,...
    AVG_EXPANSION,STD_EXPANSION,...
    AVG_EEE_COEFF_JET,STD_EEE_COEFF_JET,...
    AVG_ENTRAINMENT_JET,STD_ENTRAINMENT_JET,...
    AVG_EXPANSION_JET,STD_EXPANSION_JET,...
    PLUME_VOLUME,PLUME_HEIGHT,PLUME_TOP_RADIUS,...
    MORTON_ENTRAINMENT,...
    AVG_RELATIVE_DENSITY_JET_HEIGHT ] = preallocateVectors
%preallocateVectors Summary of this function goes here #TODO
%   Detailed explanation goes here #TODO
%
%   Last edit: Taryn Black, 11 July 2016

%%  Global variables that are used
    global nTIMESTEPS
    
%%  Preallocation for gas volume fraction
    surf_EPG = 0;%#TODO

%%  Preallocation for entrainment calculations
    AVG_EEE_COEFF = zeros(1,nTIMESTEPS);
    STD_EEE_COEFF = zeros(1,nTIMESTEPS);
    AVG_ENTRAINMENT  = zeros(1,nTIMESTEPS);
    STD_ENTRAINMENT  = zeros(1,nTIMESTEPS);
    AVG_EXPANSION  = zeros(1,nTIMESTEPS);
    STD_EXPANSION  = zeros(1,nTIMESTEPS);
    AVG_EEE_COEFF_JET = zeros(1,nTIMESTEPS);
    STD_EEE_COEFF_JET = zeros(1,nTIMESTEPS);
    AVG_ENTRAINMENT_JET  = zeros(1,nTIMESTEPS);
    STD_ENTRAINMENT_JET  = zeros(1,nTIMESTEPS);
    AVG_EXPANSION_JET  = zeros(1,nTIMESTEPS);
    STD_EXPANSION_JET  = zeros(1,nTIMESTEPS);
    PLUME_VOLUME = zeros(1,nTIMESTEPS);
    PLUME_HEIGHT = zeros(1,nTIMESTEPS);
    PLUME_TOP_RADIUS = zeros(1,nTIMESTEPS);
    MORTON_ENTRAINMENT = zeros(1,nTIMESTEPS);
    
%%  Preallocation for relative density
    AVG_RELATIVE_DENSITY_JET_HEIGHT = zeros(1,nTIMESTEPS);


end

