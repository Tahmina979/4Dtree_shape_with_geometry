clc;
clear;
close all;

% This code will registered all six 4D sequences in set 1
%To register two sequences in set2 and set3, we refer to run
%'Spatio_temporal_registration_between_two.m' script. You can register any number of 4D trees using this script with some addition/deletion of code segments.
% However, the registered 4D trees of our three different sets are saved in ‘Registered_data/Set1_reg_to1(p=0.5)’ (Set 1), ‘Registered_data/ all_lamda=1(L3,3-9), Set2’ (set 2), ‘Registered_data/ all_lamda=1(4(L3)-8(L4)), Set3’ (Set 3)
addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels");
%This code will take time as the set has complex structured trees
lam_m = 1; 
lam_s = 1;
lam_p = 1;

total_seq=6;
tree_num=1;

%We take all the 4D sequences that are spatialy registered within each sequence. We got these by
%runnung 'Spatio_temporal_registration_between_two.m' script

ctree_1 = load('Registered_data/All_within_sequence_registered_samples/seq1_L4_within.mat');

ctree_2 = load('Registered_data/All_within_sequence_registered_samples/seq2_within.mat');

ctree_5 = load('Registered_data/All_within_sequence_registered_samples/seq5_within_L4.mat');

ctree_6 = load('Registered_data/All_within_sequence_registered_samples/seq6_within_L4.mat');

ctree_7 = load('Registered_data/All_within_sequence_registered_samples/seq7_L4_within.mat');

ctree_10 = load('Registered_data/All_within_sequence_registered_samples/seq10_within.mat');

% Now we resample each 4D sequence to equal length (we choose the highest
% length among all which is the length of the first sample. So first sample
% doesn't require resampling)


%first 4D tree
num_within_seq1=numel(ctree_1.aligned_q_tree_1);
for i=1:num_within_seq1
 aligned_q_tree_1{i}=ctree_1.aligned_q_tree_1{1,i};
 sequence{tree_num}=aligned_q_tree_1{i};
 tree_num=tree_num+1;
end

%second 4D tree
num_within_seq2=numel(ctree_2.aligned_q_tree_2);
for i=1:num_within_seq2
 aligned_q_tree_2{i}=ctree_2.aligned_q_tree_2{1,i};
end
[Mu,eigenVector, EVals,tree_pca,used_qCompTrees_Ready,qX]=transfer_to_PCA_space_resamp(aligned_q_tree_2,num_within_seq2);
X1=qX(1:num_within_seq2,:)'; 
X1=ReSampleCurve(X1,7); %seven is the highest length of all curves, so we interpolate all into seven
for j=1:size(X1,2)
    aligned_q_tree_2{j}=unflattenCompTree_4layers_rad(X1(:,j)', lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
    sequence{tree_num}=aligned_q_tree_2{j};
    tree_num=tree_num+1;
end

clear X1;

%third 4D tree
num_within_seq5=numel(ctree_5.aligned_q_tree_5);
for i=1:num_within_seq5
 aligned_q_tree_5{i}=ctree_5.aligned_q_tree_5{1,i};
end
[Mu,eigenVector, EVals,tree_pca,used_qCompTrees_Ready,qX]=transfer_to_PCA_space_resamp(aligned_q_tree_5,num_within_seq5);
X1=qX(1:num_within_seq5,:)'
X1=ReSampleCurve(X1,7); %seven is the highest length of all curves, so we interpolate all into seven
for j=1:size(X1,2)
    aligned_q_tree_5{j}=unflattenCompTree_4layers_rad(X1(:,j)', lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
    sequence{tree_num}=aligned_q_tree_5{j};
    tree_num=tree_num+1;
end

clear X1;

%fourth 4D tree
num_within_seq6=numel(ctree_6.aligned_q_tree_6);
for i=1:num_within_seq6
 aligned_q_tree_6{i}=ctree_6.aligned_q_tree_6{1,i};
end

[Mu,eigenVector, EVals,tree_pca,used_qCompTrees_Ready,qX]=transfer_to_PCA_space_resamp(aligned_q_tree_6,num_within_seq6);
X1=qX(1:num_within_seq6,:)'; 
X1=ReSampleCurve(X1,7); %seven is the highest length of all curves, so we interpolate all into seven
for j=1:size(X1,2)
    aligned_q_tree_6{j}=unflattenCompTree_4layers_rad(X1(:,j)', lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
    sequence{tree_num}=aligned_q_tree_6{j};
    tree_num=tree_num+1;
end
clear X1;

%Fifth 4D tree
num_within_seq7=numel(ctree_7.aligned_q_tree_7);
for i=1:num_within_seq7
 aligned_q_tree_7{i}=ctree_7.aligned_q_tree_7{1,i};
end

[Mu,eigenVector, EVals,tree_pca,used_qCompTrees_Ready,qX]=transfer_to_PCA_space_resamp(aligned_q_tree_7,num_within_seq7);
X1=qX(1:num_within_seq7,:)'; 
X1=ReSampleCurve(X1,7); %seven is the highest length of all curves, so we interpolate all into seven
for j=1:size(X1,2)
    aligned_q_tree_7{j}=unflattenCompTree_4layers_rad(X1(:,j)', lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
    sequence{tree_num}=aligned_q_tree_7{j};
    tree_num=tree_num+1;
end

clear X1;

%Sixth 4D tree
num_within_seq10=numel(ctree_10.aligned_q_tree_10);
for i=1:num_within_seq10
 aligned_q_tree_10{i}=ctree_10.aligned_q_tree_10{1,i};
end
[Mu,eigenVector, EVals,tree_pca,used_qCompTrees_Ready,qX]=transfer_to_PCA_space_resamp(aligned_q_tree_10,num_within_seq10);
X1=qX(1:num_within_seq10,:)'; 
X1=ReSampleCurve(X1,7); %seven is the highest length of all curves, so we interpolate all into seven
for j=1:size(X1,2)
    aligned_q_tree_10{j}=unflattenCompTree_4layers_rad(X1(:,j)', lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
    sequence{tree_num}=aligned_q_tree_10{j};
    tree_num=tree_num+1;
end

%Add as  many 4D tree as you need


% Spatial registration accross all sequences
lam_m = 1; 
lam_s = 1;
lam_p = 0.5;

count=0;
for i=1:7 % total sample in a sequence which is 7 in this set(the highest length of the curves)
    control=i;
    Q2=sequence{i};
    for j=1:total_seq-1 % In set 1 total number of 4D trees are six
        control=control+7; %total sample in a sequence which is 7 in this set(the highest length of the curves)
        Q1=sequence{control};
        [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
        sequence{control}=Q1p;
        count=count+1
    end
end


num_within_seq1=7;
num_within_seq2=7;
num_within_seq5=7;
num_within_seq6=7;
num_within_seq7=7;
num_within_seq10=7;

toc

%Restoring all the curves
for i=1:num_within_seq1
  aligned_q_tree_1{i}=sequence{i};
end

p=num_within_seq1+1; p_end=num_within_seq1+num_within_seq2;

for i=p:p_end
  aligned_q_tree_2{i-p+1}=sequence{i};
end


p=p_end+1; p_end=p_end+num_within_seq5;
for i=p:p_end
   aligned_q_tree_5{i-p+1}=sequence{i};
end

p=p_end+1; p_end=p_end+num_within_seq6;
for i=p:p_end
 aligned_q_tree_6{i-p+1}=sequence{i};
end

p=p_end+1; p_end=p_end+num_within_seq7;

for i=p:p_end
   aligned_q_tree_7{i-p+1}=sequence{i};
end

p=p_end+1; p_end=p_end+num_within_seq10;
for i=p:p_end
   aligned_q_tree_10{i-p+1}=sequence{i};
end


%Perform temporal registration
lam_m = 1; 
lam_s = 1;
lam_p = 1;

clear X1;
clear X2;

%Mapping into PCA space, You need to look into and update this function if you add/delete 4D tree_shapes in/from the set
[Mu,eigenVector, EVals,tree_pca,used_qCompTrees_Ready,qX,lambda]=transfer_to_PCA_space(aligned_q_tree_1,aligned_q_tree_2,aligned_q_tree_5,aligned_q_tree_6,aligned_q_tree_7,aligned_q_tree_10,num_within_seq1,num_within_seq2,num_within_seq5,num_within_seq6,num_within_seq7,num_within_seq10);

tic

X1=tree_pca(1:num_within_seq1,:)'; %first sequence

X(:,1)=reshape(X1,[],1);

index=[num_within_seq1+1,num_within_seq1+num_within_seq2,num_within_seq1+num_within_seq2+1,num_within_seq1+num_within_seq2+num_within_seq5,num_within_seq1+num_within_seq2+num_within_seq5+1,num_within_seq1+num_within_seq2+num_within_seq5+num_within_seq6,num_within_seq1+num_within_seq2+num_within_seq5+num_within_seq6+1,num_within_seq1+num_within_seq2+num_within_seq5+num_within_seq6+num_within_seq7,num_within_seq1+num_within_seq2+num_within_seq5+num_within_seq6+num_within_seq7+1,num_within_seq1+num_within_seq2+num_within_seq5+num_within_seq6+num_within_seq7+num_within_seq10];


j=2;
for i=1:2:10 %(total_seq-1)*2=10
    X2=tree_pca(index(i):index(i+1),:)';
    [dist,X2n,q2n,X1,X2,q1, gamI, distbefore]=mygeod(X1, X2);
    X(:,j)=reshape(X2n,[],1); % Temporally alligned 4D trees stored here
    j=j+1;
  
end

%Inverse mapping
start=1;
ending=7; %total sample in a sequence which is 7 in this set(the highest length of the curves)
for i =1:total_seq 
    seq(:,start:ending)=reshape(X(:,i),[size(X1,1),size(X1,2)]);
    start=start+7;
    ending=ending+7;
end

start=1;
ending=size(X1,2);
for j=1:total_seq
    for i=start:ending
      
       tree1=seq(:,i)'*eigenVector';

       tree=tree1+Mu';
       original_data = inverseYeoJohnson(tree, lambda(i));
       recon{i}=original_data;
       %recon{i}=tree;

       sequence{i}=unflattenCompTree_4layers_rad(recon{i}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
      
    end 
    start=start+7;
    ending=ending+7;
    
end
toc

%Restoring all the curves
for i=1:num_within_seq1
  aligned_q_tree_1{i}=sequence{i};
end

p=num_within_seq1+1; p_end=num_within_seq1+num_within_seq2;

for i=p:p_end
  aligned_q_tree_2{i-p+1}=sequence{i};
end


p=p_end+1; p_end=p_end+num_within_seq5;
for i=p:p_end
   aligned_q_tree_5{i-p+1}=sequence{i};
end

p=p_end+1; p_end=p_end+num_within_seq6;
for i=p:p_end
 aligned_q_tree_6{i-p+1}=sequence{i};
end

p=p_end+1; p_end=p_end+num_within_seq7;

for i=p:p_end
   aligned_q_tree_7{i-p+1}=sequence{i};
end

p=p_end+1; p_end=p_end+num_within_seq10;
for i=p:p_end
   aligned_q_tree_10{i-p+1}=sequence{i};
end

%You can save these for future computation 
%To save in .m file, use save('seq2s_L4_alligned','aligned_q_tree_2') (an example for second sample)
% To visualize use, (an example for second sample)
%{
for i=1:num_within_seq2
    PC{i}=qCompTree_to_CompTree_rad_4layers(aligned_q_tree_2{j}); %This PC is the 4D tree in the original space
end
Then visualize PC, with the functions used in previous scripts (both in skeleton level or geometry)
%}