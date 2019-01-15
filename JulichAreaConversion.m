% JulichAreaConversion.m
%
% Convert Julich area names to brain regions
%
% Inputs:
%       - receptor_data_area_atlas_correlation.xlsx
%       - receptor_data.xls
%
% Output:
%       - julich_areas.txt (Ordered vector of brain region names for Julich data)




julich_rec_data = xlsread('.\data\Human\receptor_data.xls');

julich_area_conversion = xlsread('.\data\Human\receptor_data_area_atlas_correlation.xlsx');