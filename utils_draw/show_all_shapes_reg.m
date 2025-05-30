function show_all_shapes_reg(CT,count)
    gcf=figure;
    set(gcf,'visible', 'on');
    set(gcf, 'color', 'w');
    view(0, 90);
    axis equal; hold on;
    box on; hold on;

    set(gca,'XColor', 'none','YColor','none','ZColor','none');

% c = distinguishable_colors(70);


default_c = get(gca,'colororder');

for i=1:numel(CT)
    fprintf('Drawing %d-th/%d tree...', i, length(CT));
    tm = tic;
    spacing=2; %This can vary to adjust the spacing between two sample in a 4D sequence
    x_added = (mod(i-1, 10)) * spacing;
    if i==1
        x_added = 0;
    end

    %%Uncomment it if different spacing is needed between samples
    %{
    elseif i==2
    x_added = (mod(i-1, 10)) * 1.3;
    elseif i==3
    x_added = (mod(i-1, 10)) *1.2;
    elseif i==4
    x_added = (mod(i-1, 10)) * 1.35;
    end
    %elseif i==5
    %x_added = (mod(i-1, 10)) * 112;
    %elseif i==7
    %x_added = (mod(i-1, 10)) * 106;
    %x_added = (mod(i-1, 10)) * 126;
    %elseif i==8
    %x_added = (mod(i-1, 10)) * 122;
    %x_added = (mod(i-1, 10)) * 140;
    %elseif i==9
    %x_added = (mod(i-1, 10)) * 138;
    %x_added = (mod(i-1, 10)) * 153;
    %elseif i==10
    %x_added = (mod(i-1, 10)) * 155;
    %x_added = (mod(i-1, 10)) * 166;
    %else
        %x_added = (mod(i-1, 10)) * 80;
    %end
    %}
   y_added = 0;
  
   if count==1
        z_added=0;
   end


   % CT{i}
    clear first_pt_x first_pt_y first_pt_z;
    
    first_pt_x = (CT{i}.beta0(1, 1));
    first_pt_y = (CT{i}.beta0(2, 1));
    first_pt_z = (CT{i}.beta0(3, 1));
    
    X0 = CT{i}.beta0(1, :)+ x_added-first_pt_x ;
    Y0 = CT{i}.beta0(2, :) + y_added-first_pt_y ;
    Z0 = CT{i}.beta0(3, :)+z_added- first_pt_z;
    %R0 = CT{i}.beta0_rad;
   
    plot3(X0, Y0, Z0, 'Color',default_c(2,:), 'LineWidth', 1);
 
    for j=1: numel(CT{i}.beta)
        clear X Y Z
        
        X = CT{i}.beta{j}(1, :)+ x_added- first_pt_x ;
        Y = CT{i}.beta{j}(2, :) + y_added- first_pt_y ;
        Z = CT{i}.beta{j}(3, :)+z_added- first_pt_z;
        %R= CT{i}.beta_rad{j};
        
        plot3(X, Y, Z, 'Color',default_c(1,:), 'LineWidth', 1); hold on;  
      
        for k = 1: numel(CT{i}.beta_children{j}.beta)
            clear X Y Z
            X = CT{i}.beta_children{j}.beta{k}(1, :)+ x_added- first_pt_x ;
            Y = CT{i}.beta_children{j}.beta{k}(2, :)+ y_added- first_pt_y ;
            Z = CT{i}.beta_children{j}.beta{k}(3, :)+z_added- first_pt_z;
            %R= CT{i}.beta_children{j}.beta_rad{k};
            
            plot3(X, Y, Z, 'Color',default_c(3,:), 'LineWidth', 1); hold on;
            
            for t = 1: numel(CT{i}.beta_children{j}.beta_children{k}.beta)
                clear X Y Z
                X = CT{i}.beta_children{j}.beta_children{k}.beta{t}(1, :)+ x_added- first_pt_x ;
                Y = CT{i}.beta_children{j}.beta_children{k}.beta{t}(2, :) + y_added- first_pt_y ;
                Z = CT{i}.beta_children{j}.beta_children{k}.beta{t}(3, :)+z_added- first_pt_z;
                % R = CT{i}.beta_children{j}.beta_children{k}.beta_rad{t};
                plot3(X, Y, Z, 'Color',default_c(4,:), 'LineWidth', 1); hold on;
            
            end
            
        end
        
    end
    
    T_cost = toc(tm);
    fprintf('-Done, time cost: %.2f\n', T_cost)
        
end

end