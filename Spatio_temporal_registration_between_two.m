
clear;
clc;
close all;

addpath('utils_data','utils_draw','utils_funcs',"utils_statModels","OpenCurvesRn");
%These parameters value range from 0 to 1 and control spatial registration. Default: all to 1

lam_m = 1; 
lam_s = 1;
lam_p = 1;

data_path5= 'NeuroData/Sample5/'; %Can change the sample
[all_qCompTrees5, all_compTrees5] = load_botanTrees_rad(data_path5);

show_all_shapes_sp_reg_l4(all_compTrees5,1);


num_within_seq5=numel(all_qCompTrees5);
%Within sequence registration
aligned_q_tree_5{1}=all_qCompTrees5{1};
for i=2:num_within_seq5
    
    Q1=all_qCompTrees5{i};
    Q2 = aligned_q_tree_5{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_5{i}=Q1p;
end
%save("seq5_within.mat", "aligned_q_tree_5"); %can save intermediate results

data_path6= 'NeuroData/Sample6/'; %Can chane the sample
[all_qCompTrees6, all_compTrees6] = load_botanTrees_rad(data_path6);

show_all_shapes_sp_reg_l4(all_compTrees6,2);

num_within_seq6=numel(all_qCompTrees6);
%Within sequence registration
aligned_q_tree_6{1}=all_qCompTrees6{1};
for i=2:num_within_seq6
    
    Q1=all_qCompTrees6{i};
    Q2 = aligned_q_tree_6{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_6{i}=Q1p;
end
%save("seq6_within.mat", "aligned_q_tree_6"); %can save intermediate results

% You can add as many samples as you like. We performed within
% sequence registration for six samples in set 1 and two samples for set2
% and set3. All the results are saved in
% Registered_Data/All_within_sequence_registered_samples folder


%{
data_path7= 'NeuroData/Sample7/';
[all_qCompTrees7, all_compTrees7] = load_botanTrees_rad(data_path7);
num_within_seq7=numel(all_qCompTrees7);

aligned_q_tree_7{1}=all_qCompTrees7{1};
for i=2:num_within_seq7
    count=count+1
    Q1=all_qCompTrees7{i};
    Q2 = aligned_q_tree_7{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_7{i}=Q1p;
end
%save("seq7_within.mat", "aligned_q_tree_7"); %can save intermediate results


data_path4= 'NeuroData/Sample8/';
[all_qCompTrees8, all_compTrees8] = load_botanTrees_rad(data_path4);
num_within_seq8=numel(all_qCompTrees8);

aligned_q_tree_8{1}=all_qCompTrees8{1};
for i=2:num_within_seq8
    count=count+1
    Q1=all_qCompTrees8{i};
    Q2 = aligned_q_tree_8{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_8{i}=Q1p;
end
save("seq8_within.mat", "aligned_q_tree_8");

data_path5= 'NeuroData/Sample9/';
[all_qCompTrees9, all_compTrees9] = load_botanTrees_rad(data_path5);
num_within_seq9=numel(all_qCompTrees9);

aligned_q_tree_9{1}=all_qCompTrees9{1};
for i=2:num_within_seq9
    count=count+1
    Q1=all_qCompTrees9{i};
    Q2 = aligned_q_tree_9{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_9{i}=Q1p;
end
save("seq9_within.mat", "aligned_q_tree_9");


data_path6= 'NeuroData/Sample10/';
[all_qCompTrees10, all_compTrees10] = load_botanTrees_rad(data_path6);
num_within_seq10=numel(all_qCompTrees10);

aligned_q_tree_10{1}=all_qCompTrees10{1};
for i=2:num_within_seq10
    count=count+1
    Q1=all_qCompTrees10{i};
    Q2 = aligned_q_tree_10{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_10{i}=Q1p;
end
save("seq10_within.mat", "aligned_q_tree_10");
%}

%Resample 4D sequences to equal length (We choose the highest length in two/all sequences)
[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready,qX,lambda]=transfer_to_PCA_geod_space(aligned_q_tree_5,aligned_q_tree_6,num_within_seq5,num_within_seq6); % for two sequences

X1=qX(1:num_within_seq5,:)'; %first sequence
X2=qX(num_within_seq5+1:num_within_seq5+num_within_seq6,:)'; % second sequence

%For 4D sequence/curve handling set these parameters in default values
lam_m = 1; 
lam_s = 1;
lam_p = 1;

if num_within_seq5>num_within_seq6
    number_of_seq=num_within_seq5;
    X2=ReSampleCurve(X2,num_within_seq5); %second one resampled
    for j=1:size(X2,2)
         aligned_q_tree_6{j}=unflattenCompTree_4layers_rad(X2(:,j)', lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
    end
        num_within_seq6=num_within_seq5;
elseif num_within_seq6>num_within_seq5
    
    number_of_seq=num_within_seq6;
    X1=ReSampleCurve(X1,num_within_seq6); % first one resampled
    for j=1:size(X1,2)
         aligned_q_tree_5{j}=unflattenCompTree_4layers_rad(X1(:,j)', lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
    end
        num_within_seq5=num_within_seq6;
end




tic
count=0
%spatial registration across_sequence for two sequences
for i=1:number_of_seq
    count=count+1
    Q1=aligned_q_tree_6{i};
    Q2=aligned_q_tree_5{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2_within(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_6{i}=Q1p;
    aligned_q_tree_5{i}=Q2p
end

toc

% Perform temporal registration between two

[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready,qX,lambda]=transfer_to_PCA_geod_space(aligned_q_tree_5,aligned_q_tree_6,num_within_seq5,num_within_seq6);

X1=tree_pca(1:num_within_seq5,:)';
X2=tree_pca(num_within_seq5+1:num_within_seq5+num_within_seq6,:)';


[dist,X2n,X1,X2,gamI, distbefore]=mygeod(X1,X2); %Perfrom the temporal alignment, X2n is the temporally aligned 4D sequence/curve of X2(second sample) onto X1(first sample) (in PCA space)

%inverse PCA mapping
for i=1:size(X2n,2)
       tree1=X2n(:,i)'*eigenVector';
       tree=tree1+Mu';
       original_data = inverseYeoJohnson(tree, lambda(num_within_seq5+i));
       recon{i}= original_data;
end

for j=1:size(X2n,2)
  
   aligned_q_tree_6{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{j});
   PC{j}=qCompTree_to_CompTree_rad_4layers(aligned_q_tree_5{j}); %for visualization
   PC1{j}=qCompTree_to_CompTree_rad_4layers(aligned_q_tree_6{j}); %for visualization

end


show_all_shapes_sp_reg_l4(PC,3);
show_all_shapes_sp_reg_l4(PC1,4);

%% We can save the results  for future computation, We save for some pairs in Registerd_data folder
%save("seq6_reg.mat", "aligned_q_tree_6");
%save("seq5_reg.mat", "aligned_q_tree_5");
