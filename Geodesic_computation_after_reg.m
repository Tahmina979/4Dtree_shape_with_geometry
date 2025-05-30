clc;
clear;
close all;

addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels");

lam_m = 1; 
lam_s = 1;
lam_p = 1;


% we compute geodesic by loading data from saved registered samples.
% Sequences should be spatiotemporally registered in step2.

%Load from saved registered samples (first sample)
ctree_1 = load('Registered_data/all_lamda=1(L4,5,10)/seq5_reg.mat');

num_within_seq1=numel(ctree_1.aligned_q_tree_5); %name 'aligned_q_tree_1' comes from mat file 

for i=1:num_within_seq1
    q_tree_1{i}=ctree_1.aligned_q_tree_5{i};
end

%Load from saved registered samples (second sample)
ctree_2 = load('Registered_data/all_lamda=1(L4,5,10)/seq10_reg.mat');
num_within_seq2=numel(ctree_2.aligned_q_tree_10);

for i=1:num_within_seq2
    q_tree_2{i}=ctree_2.aligned_q_tree_10{i};
end

[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready,qX,lambda]=transfer_to_PCA_geod_space(q_tree_1,q_tree_2,num_within_seq1,num_within_seq2);

X1=tree_pca(1:num_within_seq1,:)';
X2=tree_pca(num_within_seq1+1:num_within_seq1+num_within_seq2,:)';

%geod computation

tic
Xgeod = computeGeodesic(X1,X2);
toc
count=0;

nSamples=7; %How many intermediate samples we want to take on the straightline between two 4D trees
lambda1=lambda(1:num_within_seq1);
lambda2=lambda(num_within_seq1+1:num_within_seq1+num_within_seq2);

for k=1:size(Xgeod,1)
  t = (k-1.0)/(nSamples-1);
  in_lambda=(1-t)*lambda1 + t*lambda2;
   for i=1:size(X2,2)
       tree1=Xgeod(k,:,i)*eigenVector';
       tree=tree1+Mu';
       original_data = inverseYeoJohnson(tree, in_lambda(i));
       recon{i}=original_data;
       %recon{i}=tree;
   end
   
for j=1:size(X2,2)
   PC1{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
   all_compTrees1{j}=qCompTree_to_CompTree_rad_4layers(PC1{j});

end
   
  show_all_shapes_reg(all_compTrees1,1);
%%uncomment the code segment if you want to see the rendered geometry of
%%the tree models (rendering will take 5-10 minutes)


    idx1=k;
    idx2=k;
    data_path='check'%'geod_between_4Dsamples_after_Reg'
    % % ===== Save Objs =====
    addpath('GetOBJ')
    obj_folder = saveGeoObjs_compTrees_rad_4layers(all_compTrees1, data_path, idx1, idx2);
    
    % ===== Render Objs =====
    addpath('RenderOBJ')
    addpath('RenderOBJ/func_render/');
    addpath('RenderOBJ/func_other/');
    run renderSequence_geod.m;


end

