clear
close all


% setting up the initial configuration 
P.y0 = 180;%nm, corresponds with spectrin bundle length
P.x0 = sqrt(P.y0^2 - (P.y0/2)^2);%nm; height of the triangles
P.max_x0 = P.x0*15;
P.max_y0 = (ceil(P.max_x0/P.y0)-1)*P.y0;

% create triangular mesh
[r_s,T_s,edges_s,edge_type] = create_mesh(P);
% r_s:position of the nodes, edges_s: edge connectivity matrix, edge_type :
% type of edge
P.actin = unique(edges_s);%nodes that correspond to short F-actin
P.rope_length = 360;%initial length of connecting edges

% attach connecting edges to the mesh
% P.rope correspond to the focal adhesion nodes
% P.link correspond to link nodes
aux_n = size(r_s,1);
aux = find(r_s(P.actin,1) == max(r_s(P.actin,1)));
r_s = [r_s; (r_s(P.actin(aux),1)+P.rope_length) r_s(P.actin(aux),2)];
edges_s = [edges_s; P.actin(aux) aux_n+(1:length(aux))'];
edge_type = [edge_type;ones(size(aux))];
P.rope = (aux_n+1:size(r_s,1))';
P.link = P.actin(aux);
P.actin(aux) = [];

aux_n = size(r_s,1);
aux = find(r_s(P.actin,1) == min(r_s(P.actin,1)));
r_s = [r_s; (r_s(P.actin(aux),1)-P.rope_length) r_s(P.actin(aux),2)];
edges_s = [edges_s; P.actin(aux) aux_n+(1:length(aux))'];
edge_type = [edge_type;ones(size(aux))];
P.rope = [P.rope;(aux_n+1:size(r_s,1))'];
P.link = [P.link;P.actin(aux)];
P.actin(aux) = [];

aux_n = size(r_s,1);
aux = find(r_s(P.actin,2) == max(r_s(P.actin,2)) & ...
    (r_s(P.actin,1) < max(r_s(P.actin,1)) & r_s(P.actin,1) > min(r_s(P.actin,1))));
r_s = [r_s; r_s(P.actin(aux),1) (r_s(P.actin(aux),2)+P.rope_length)];
edges_s = [edges_s; P.actin(aux) aux_n+(1:length(aux))'];
edge_type = [edge_type;ones(size(aux))];
P.rope = [P.rope;(aux_n+1:size(r_s,1))'];
P.link = [P.link;P.actin(aux)];
P.actin(aux) = [];

%  
aux_n = size(r_s,1);
aux = find(r_s(P.actin,2) == min(r_s(P.actin,2)) & ...
    (r_s(P.actin,1) < max(r_s(P.actin,1)) & r_s(P.actin,1) > min(r_s(P.actin,1))));
r_s = [r_s; r_s(P.actin(aux),1) (r_s(P.actin(aux),2)-P.rope_length)];
edges_s = [edges_s; P.actin(aux) aux_n+(1:length(aux))'];
edge_type = [edge_type;ones(size(aux))];
P.rope = [P.rope;(aux_n+1:size(r_s,1))'];
P.link = [P.link;P.actin(aux)];
P.actin(aux) = [];

% change patch to 3D
r_s = [r_s zeros(size(r_s,1),1)];
r_s(P.rope,3) = -10*ones(size(P.rope));


% % % % create a plot to check if the mesh is correct
% figure;
% trimesh(T_s,r_s(:,1),r_s(:,2),r_s(:,3),'edgecolor',[0.75 0.75 0.75])
% hold on 
% % plot spectrin bundles
% aux_e = find(edge_type == 0);
% for l=1:size(aux_e,1)
%     plot3([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
%        [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
%        [r_s(edges_s(aux_e(l),1),3),r_s(edges_s(aux_e(l),2),3)],'b')
% end
% aux_e = find(edge_type == 1);
% for l=1:size(aux_e,1)
%     plot3([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
%        [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
%        [r_s(edges_s(aux_e(l),1),3),r_s(edges_s(aux_e(l),2),3)],'k')
% end
% axis('equal')
% plot3(r_s(P.actin,1),r_s(P.actin,2),r_s(P.actin,3),'ms')%plot actin short filaments 
% plot3(r_s(P.rope,1),r_s(P.rope,2),r_s(P.rope,3),'ko')%
% plot3(r_s(P.link,1),r_s(P.link,2),r_s(P.link,3),'rv')%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% define parameters
P.k0 = 1;%pN/nm,spectrin spring constant
P.d00 = P.y0;%nm, spectrin resting length

P.k1 = 1;%pN/nm,rope spring constant
P.d01 = P.rope_length-P.rope_length/4;%nm, spectrin resting length

P.delta_t =0.002;% 0.002;%s, time step size
P.zeta = 1.25;%pN s/nm, drag

P.t_ini=0;%s, initial time
P.t_end = 180;%s,end time
aux_t = P.t_ini:P.delta_t:P.t_end;

% save variables
save_aux = 10;%saving frequency
r_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);

ll=1;
r_save{ll,1} = r_s;

for k=2:length(aux_t)
%     force vector
    f = zeros(size(r_s));

   
    for l = 1:size(edges_s,1) 
%         calculate the distance of edges
        r_ij = r_s(edges_s(l,1),:) - r_s(edges_s(l,2),:);
        d = sqrt(dot(r_ij,r_ij,2));
%         depending of the type of edge compute the force and store the
%         value

%         for spring edges
        if edge_type(l) == 0
            f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.k0*r_ij*(d-P.d00)/d;
            f(edges_s(l,2),:) =  f(edges_s(l,2),:) + P.k0*r_ij*(d-P.d00)/d;
%         for connecting edges
        elseif edge_type(l) == 1
            f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.k1*r_ij*(d-P.d01)/d;
            f(edges_s(l,2),:) =  f(edges_s(l,2),:) + P.k1*r_ij*(d-P.d01)/d;
        end
    end
%     update the nodes position 
    r_s(P.actin,:) = P.delta_t*f(P.actin,:)./P.zeta + r_s(P.actin,:);
    r_s(P.link,:) = P.delta_t*f(P.link,:)./P.zeta + r_s(P.link,:);
  
%     save data 
    if mod(k,save_aux) == 1
        ll=ll+1;
        r_save{ll,1} = r_s;       
    end
end
