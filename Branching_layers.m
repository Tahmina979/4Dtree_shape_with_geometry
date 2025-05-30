clc;
clear;
close all;

addpath('utils_draw', 'utils_data');


data_path= 'NeuroData/Sample1_L4/'; %Change with the sample you want to visualize

[all_qCompTrees1, all_compTrees1] = load_botanTrees_rad(data_path);

show_all_shapes_reg(all_compTrees1,1);
