clc;
clear;
close all;

addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels",'Registered_tomato_plants');

lam_m = 1; 
lam_s = 1;
lam_p = 1;
Number_of_sequences=6;
Number_of_sample=7; % number of sample in a sequence


%We load spatiotemporally registered 4D tree models from set 1

ctree_1 = load('Registered_data/Set1_reg_to1(p=0.5)/seq1_L4_within.mat');
ctree_2 = load('Registered_data/Set1_reg_to1(p=0.5)/seq2s_L4_alligned.mat');
ctree_5 = load('Registered_data/Set1_reg_to1(p=0.5)/seq5s_L4_alligned.mat');
ctree_6 = load('Registered_data/Set1_reg_to1(p=0.5)/seq6s_L4_alligned.mat');
ctree_7 = load('Registered_data/Set1_reg_to1(p=0.5)/seq7s_L4_alligned.mat');
ctree_10 = load('Registered_data/Set1_reg_to1(p=0.5)/seq10s_L4_alligned.mat');

num_within_seq1=numel(ctree_1.aligned_q_tree_1);
num_within_seq2=numel(ctree_2.aligned_q_tree_2);
num_within_seq5=numel(ctree_5.aligned_q_tree_5);
num_within_seq6=numel(ctree_6.aligned_q_tree_6);
num_within_seq7=numel(ctree_7.aligned_q_tree_7);
num_within_seq10=numel(ctree_10.aligned_q_tree_10);


[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready,qX,lambda]=transfer_to_PCA_space_mean_mode(ctree_1,ctree_2,ctree_5,ctree_6,ctree_7,ctree_10,num_within_seq1,num_within_seq2,num_within_seq5,num_within_seq6,num_within_seq7,num_within_seq10);


%reshaping curves
start=1;
for i=1:Number_of_sequences 
    ending=start+Number_of_sample-1; 
    X(:,i)=reshape(tree_pca(start:ending,:)',[],1);
    start=ending+1;
end

control=1;
% Mu_s is the mean curve
for i=1:Number_of_sample
    sum=0;
    for j=i:Number_of_sample:Number_of_sample*Number_of_sequences  
        sum=sum+lambda(j);
    end
    mean_lambda(i)=sum/Number_of_sample;
    
end

%PCA on curves
[Mu_s, eigenVectors_s, EVals_s] = performEigenAnalysis(X);

alpha_1 = sqrt(EVals_s(1));
alpha_2 = sqrt(EVals_s(2));
alpha_3 = sqrt(EVals_s(3));
alpha_4 = sqrt(EVals_s(4));

Range=1.5; Step=0.5;

std=[-1.5,-1.0,-0.5,0,0.5,1.0,1.5];

alpha_1_vec = -Range*alpha_1: (Step*alpha_2): Range*alpha_1; %Take alpha_2 to see the variations in the second princiapal direction and so on.

mode_length= numel(alpha_1_vec);

tic
for i=1:mode_length
    qX_PC1{i} = alpha_1_vec(i)*eigenVectors_s(:,1) + Mu_s;
end
toc

X=zeros(mode_length,size(tree_pca,2),Number_of_sample); 

for i=1:mode_length
    X(i,:,:)=reshape(qX_PC1{i},[size(tree_pca,2),Number_of_sample]);
end

for md=1:mode_length
    mode(:,1:Number_of_sample)=X(md,:,:);
    for i=1:size(X,3)
           tree1=X(md,:,i)*eigenVector'; 
           tree=tree1+Mu';
           original_data = inverseYeoJohnson(tree, mean_lambda(i));
           recon{i}=original_data;
           %recon{i}=tree;
    end

for j=1:size(recon,2)
   mean_mode{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
  all_compTrees1{j}=qCompTree_to_CompTree_rad_4layers(mean_mode{j});
end

show_all_shapes_reg(all_compTrees1,1); %uncomment if you want to see at the skeleton level
%Rendering will take time. If want to see the geometry, uncommment the code
%segment below
%{
    data_path='check_modes_mean';

    idx1=md;
    idx2=md;
       % % ===== Save Objs =====
    addpath('GetOBJ') 
    obj_folder = saveGeoObjs_compTrees_rad_4layers(all_compTrees1, data_path, idx1, idx2);
    
    % ===== Render Objs =====
    addpath('RenderOBJ')
    addpath('RenderOBJ/func_render/');
    addpath('RenderOBJ/func_other/');
    run renderSequence_mean_mode.m
%}

end

return;
