function [ BULK_DENSITY,RELATIVE_DENSITY,AVG_RELATIVE_DENSITY_JETHEIGHT ]...
    = calculateDensities( EPG,EPS1,EPS2,EPS3,RO_G,JET_HEIGHT )
%calculateDensities Summary of this function goes here #TODO
%   Detailed explanation goes here #TODO
%
%   Last edit: Taryn Black, 21 July 2016

%%  Global variables that are used
    global RO_S1 RO_S2 RO_S3
    global PLUME_EDGE
    global IMAX JMAX KMAX GHOSTCELLS_IJK

%%  Calculate bulk density of flow
    BULK_DENSITY = (EPG.*RO_G)+(EPS1*RO_S1)+(EPS2*RO_S2)+(EPS3*RO_S3);
    
%%  Calculate density of flow relative to atmosphere
 %  Determine whether domain elements are inside or outside of plume
    IN_FLOW = EPG <= PLUME_EDGE;
    IN_ATMOS = EPG > PLUME_EDGE;
    
 %  Separate bulk density matrix into flow and atmospheric components
    FLOW_DENSITY = BULK_DENSITY.*IN_FLOW;
    ATMOS_DENSITY = BULK_DENSITY.*IN_ATMOS;
    
 %  Calculate average atmospheric density profile at each altitude gridpoint
    AVG_ATMOS_DENSITY_PROFILE = zeros(1,JMAX - GHOSTCELLS_IJK(2));
    for k = 1:(JMAX - GHOSTCELLS_IJK(2))
        ATMOS_DENSITY_AT_Z = ATMOS_DENSITY(:,:,k);
        AVG_ATMOS_DENSITY_PROFILE(k) = mean(ATMOS_DENSITY_AT_Z(ATMOS_DENSITY_AT_Z~=0));
    end
    
 %  Calculate relative density (of flow to atmosphere)
 %  Positive: flow is less dense than atmosphere (buoyant)
 %  Negative: flow is denser than atmosphere (not buoyant)
    AVG_ATMOS_DENSITY_PROFILE_3D = repmat(AVG_ATMOS_DENSITY_PROFILE',1,...
        KMAX - GHOSTCELLS_IJK(3),IMAX - GHOSTCELLS_IJK(1));
    AVG_ATMOS_DENSITY_PROFILE_3D = permute(AVG_ATMOS_DENSITY_PROFILE_3D,[3 2 1]);
    AVG_ATMOS_DENSITY_PROFILE_3D_AT_FLOW = AVG_ATMOS_DENSITY_PROFILE_3D.*IN_FLOW;
    RELATIVE_DENSITY = AVG_ATMOS_DENSITY_PROFILE_3D_AT_FLOW - FLOW_DENSITY;
    
 %  Calculate average relative density of flow at jet height
    FLOW_AT_JETHEIGHT = IN_FLOW(:,:,round(JET_HEIGHT));
    RELATIVE_DENSITY_AT_JETHEIGHT = RELATIVE_DENSITY(:,:,round(JET_HEIGHT));
    AVG_RELATIVE_DENSITY_JETHEIGHT = mean(RELATIVE_DENSITY_AT_JETHEIGHT(FLOW_AT_JETHEIGHT));


end

