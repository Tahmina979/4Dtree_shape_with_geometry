clc;
clear;
close all;

addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels",'Registered_maize_plants');

lam_m = 1; 
lam_s = 1;
lam_p = 1;
Total_seq=2;

% we load spatiotemporally registered 4D sequences from the saved one

% this script is for getting mean and modes of 4D trees in set 3. To get mean and mode from set 2, replace the instructed lines
number_of_samples=5; % set number_of_samples=6 for set 2

ctree_4 = load('Registered_data/all_lamda=1(4(L3)-8(L4)),Set3/seq4_reg.mat'); %replace with 'Registered_data/all_lamda=1(L3,3-9),Set2/seq3_reg.mat'
num_within_seq4=numel(ctree_4.aligned_q_tree_4); %replace 'aligned_q_tree_4' with 'aligned_q_tree_3'
for i=1:num_within_seq4 
    aligned_q_tree_4{i}=ctree_4.aligned_q_tree_4{1,i}; %replace 'aligned_q_tree_4' with 'aligned_q_tree_3'
    
end

ctree_8 = load('Registered_data/all_lamda=1(4(L3)-8(L4)),Set3/seq8_reg.mat'); %replace with 'Registered_data/all_lamda=1(L3,3-9),Set2/seq9_reg.mat'
num_within_seq8=numel(ctree_8.aligned_q_tree_8); %replace 'aligned_q_tree_8' with 'aligned_q_tree_9'
for i=1:num_within_seq8
    aligned_q_tree_8{i}=ctree_8.aligned_q_tree_8{1,i}; %replace 'aligned_q_tree_8' with 'aligned_q_tree_9'    
end


[Mu,eigenVector, EVals,tree_pca,used_qCompTrees_Ready,qX,lambda]=transfer_to_PCA_geod_space(aligned_q_tree_4,aligned_q_tree_8,num_within_seq4,num_within_seq8);
%reshaping curves
start=1;
for i=1:Total_seq 
    ending=start+number_of_samples-1; 
    X(:,i)=reshape(tree_pca(start:ending,:)',[],1);
    start=ending+1;
end
control=1;
% Mu_s is the mean curve
for i=1:number_of_samples
    sum=0;
    for j=i:number_of_samples:number_of_samples*Total_seq
        sum=sum+lambda(j);
    end
    mean_lambda(i)=sum/Total_seq;
    
end
%% compute modes of variation 
%PCA on curves

[Mu_s, eigenVectors_s, EVals_s] = performEigenAnalysis(X);


alpha_1 = sqrt(EVals_s(1));

Range=1.5; Step=0.5;

alpha_1_vec = -Range*alpha_1: (Step*alpha_1): Range*alpha_1;

mode_length= numel(alpha_1_vec);
std=[-1.5,-1.0,-0.5,0,0.5,1.0,1.5];
tic
for i=1:mode_length
    qX_PC1{i} = alpha_1_vec(i)*eigenVectors_s(:,1) + Mu_s;
end
toc

X=zeros(mode_length,size(tree_pca,2),num_within_seq4); 

for i=1:mode_length
    X(i,:,:)=reshape(qX_PC1{i},[size(tree_pca,2),num_within_seq4]);
end

for md=1:mode_length
    mode(:,1:number_of_samples)=X(md,:,:);
    
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

end
return;
