clc; 
%clear; 
%close all;

mesh_dir = [obj_folder, '/'];

% load the shape
N=numel(all_compTrees1);
S = cell(1, N);
for i = 1:N
    S{i} = MESH_IO.read_shape([mesh_dir, num2str(i),'.obj']);
end

%% Draw geodesic
 renderOptions = {'RotationOps',{[-90,0,0],[0,0,0]},...  
                  'CameraPos',[-5, 10],...             % default: [-10, 10]
                  'FaceAlpha',0.9,...
                  'BackgroundColor',[0.9, 0.9, 0.9]}; % you can change the background color here
              

   
    gcf=figure;
    axis equal; hold on;
    set(gcf,'visible', 'on');
    set(gcf, 'color', 'w');
    default_c = get(gca,'colororder');
   
    set(gca,'XColor', 'none','YColor','none','ZColor','none');
  

for i=1:N
    

    M = S{i};
    color_botanTree = [0.2, 0.1, 0];
    M_col = repmat(color_botanTree, M.nv, 1);    

    moveVec = [(i-1)*2.5, 0, 0];  % Spacing can vary for different 4D sequences
  
    if i==1
        moveVec = [0,0, 0];
    end
    %This segment can be used for different spacing between two consecutive samples.
    %{
    if i==4
        %moveVec = [(i-1)*1.65,0,0];
    end
   
    if  i==5
        %moveVec = [(i-1)*1.85,0, 0];
    end
    
    if i==4
        moveVec = [(i-1)*1,-1, 0];
    end
    if i==5
        moveVec = [(i-1)*1,-1.5, 0];
    end
    if i==6
        moveVec = [(i-1)*1.25, -1.6, 0];
    end
    if i==7
        moveVec = [(i-1)*1.33,-1.5, 0];
    end
     %}
    tm1 = tic();
    [~,~, S1_new] = render_mesh_sequence(M,'MeshVtxColor',M_col,...
                               'VtxPos',M.surface.VERT + repmat(moveVec, M.nv,1),... % translate the second shape such that S1 and S2 are not overlapped
                                renderOptions{:});
    hold on;

    % --- Draw all edges ---
    all_edges = get_edge_list(S{i});  % get the edge list
    vertex = S1_new.surface.VERT';    % get the vertex positions after the rotations!

    
    t = trimesh(S1_new.surface.TRIV, vertex(1,:), vertex(2,:), vertex(3,:),...
                'FaceColor',[0.5, 0.5, 0.5],...
                'EdgeColor', color_botanTree,...
                'LineWidth', 0.1);

    hold on;
    T_render= toc(tm1);
    fprintf('%d-th tree rendering done, time: %f\n', i, T_render);
    
end
