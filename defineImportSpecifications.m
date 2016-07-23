function [ IMPORTSPEC_EPG,IMPORTSPEC_EPS1,IMPORTSPEC_EPS2,IMPORTSPEC_EPS3,...
    IMPORTSPEC_UG,IMPORTSPEC_VG,IMPORTSPEC_WG,...
    IMPORTSPEC_TG,IMPORTSPEC_ROG,...
    IMPORTSPEC_VS1,IMPORTSPEC_VS2,IMPORTSPEC_VS3] = defineImportSpecifications
%defineImportSpecifications Summary of this function goes here #TODO
%   Specifies columns to read (%f) or skip (%*f) for each data file.
%
%   Last edit: Taryn Black, 11 July 2016

%%  Gas and particle volume fractions
%   File: EP_t**.txt
%   File structure:
%   1       2       3       4       5       6       7
%   EP_G    EP_S1   EP_S2   EP_S3   X       Y       Z
    IMPORTSPEC_EPG  = '%f%*f%*f%*f%*f%*f%*f';
    IMPORTSPEC_EPS1 = '%*f%f%*f%*f%*f%*f%*f';
    IMPORTSPEC_EPS2 = '%*f%*f%f%*f%*f%*f%*f';
    IMPORTSPEC_EPS3 = '%*f%*f%*f%f%*f%*f%*f';
    
%%  Directional components of gas velocity
%   File: U_G_t**.txt
%   File structure:
%   1       2       3       4       5       6
%   U_G     V_G     W_G     X       Y       Z
    IMPORTSPEC_UG = '%f%*f%*f%*f%*f%*f';
    IMPORTSPEC_VG = '%*f%f%*f%*f%*f%*f';
    IMPORTSPEC_WG = '%*f%*f%f%*f%*f%*f';
    
%%  Gas temperature
%   File: T_G_t**.txt
%   File structure:
%   1       2       3       4
%   T_G     X       Y       Z
    IMPORTSPEC_TG = '%f%*f%*f%*f';
    
%%  Gas density
%   File: Current_Density_t**.txt
%   File structure:
%   1      2      3                         4                   5   6   7
%   RO_C   RO_G   RO_C/avg_RO_air(height)   avg_RO_air(height)  X   Y   Z
    IMPORTSPEC_ROG = '%*f%f%*f%*f%*f%*f%*f';
    
%%  Vertical velocities of each particle phase
    IMPORTSPEC_VS1 = '%f%*f%*f';
    IMPORTSPEC_VS2 = '%*f%f%*f';
    IMPORTSPEC_VS3 = '%*f%*f%f';


end

