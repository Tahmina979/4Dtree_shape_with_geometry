clc;
clear;
close all;

addpath('utils_data\');


data_path= 'NeuroData/Sample1_L4/'; %Change with the sample you want to visualize

[all_qCompTrees1, all_compTrees1] = load_botanTrees_rad(data_path);

 %Saving the models as mesh in 'outputedOBJs/geodOBJs_data_path_mesh-idx1andidx2'   
 idx1=1;
 idx2=2;
 data_path_mesh='sample1_l4'
 addpath('GetOBJ');
 obj_folder = saveGeoObjs_compTrees_rad_4layers(all_compTrees1, data_path_mesh, idx1, idx2);

 %Render saved meshes for visualization
 addpath('RenderOBJ')
 addpath('RenderOBJ/func_render/');
 addpath('RenderOBJ/func_other/');
 run renderSequence.m;
