clc;
clear;
close all;

addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels");

lam_m = 1; 
lam_s = 1;
lam_p = 1;


% we compute geodesic on the raw data, if unequal length, we take interpolated one to make equal
% lenth 4D tree

%First sample
data_path1='NeuroData/Sample5/';
[q_tree_1, q_ctree_1] = load_botanTrees_rad(data_path1);
num_within_seq1=numel(q_ctree_1);


%Second sample
data_path2 = 'NeuroData/Sample10/'
[q_tree_2, q_ctree_2] = load_botanTrees_rad(data_path2);
num_within_seq2=numel(q_ctree_2);

[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready,qX,lambda]=transfer_to_PCA_geod_space(q_tree_1,q_tree_2,num_within_seq1,num_within_seq2);

X1=qX(1:num_within_seq1,:)'; 
X2=qX(num_within_seq1+1:num_within_seq1+num_within_seq2,:)';

if num_within_seq1>num_within_seq2
    num_within_seq=num_within_seq1;
    X2=ReSampleCurve(X2,num_within_seq); %second one resampled
    for j=1:size(X2,2)
         q_tree_2{j}=unflattenCompTree_4layers_rad(X2(:,j)', lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
    end
        num_within_seq2=num_within_seq;

elseif num_within_seq2>num_within_seq1
    num_within_seq=num_within_seq2;
    X1=ReSampleCurve(X1,num_within_seq); % first one resampled
    for j=1:size(X1,2)
         q_tree_1{j}=unflattenCompTree_4layers_rad(X1(:,j)', lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
    end
        num_within_seq1=num_within_seq;
end


[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready,qX,lambda]=transfer_to_PCA_geod_space(q_tree_1,q_tree_2,num_within_seq1,num_within_seq2);

X1=tree_pca(1:num_within_seq1,:)';
X2=tree_pca(num_within_seq1+1:num_within_seq1+num_within_seq2,:)';

%geod computation

tic
Xgeod = computeGeodesic(X1,X2);
toc
count=0;

nSamples=7; % How many intermediate samples we want to take on the straightline between two 4D trees
in_lambda=lambda;

for k=1:size(Xgeod,1)
  t = (k-1.0)/(nSamples-1);
  
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
    data_path='check'%'geod_between_4Dsamples_before_Reg'
    % % ===== Save Objs =====
    addpath('GetOBJ')
    obj_folder = saveGeoObjs_compTrees_rad_4layers(all_compTrees1, data_path, idx1, idx2);
    
    % ===== Render Objs =====
    addpath('RenderOBJ')
    addpath('RenderOBJ/func_render/');
    addpath('RenderOBJ/func_other/');
    run renderSequence_geod.m;


end

