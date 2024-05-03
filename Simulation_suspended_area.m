clear
close all
rng(040522)

% load spherical shape 
load('sphere.mat')
r_s = 1.06*r_s;
% setting up the initial configuration 
P.y0 = 180;%nm, corresponds with spectrin bundle length
P.x0 = sqrt(P.y0^2 - (P.y0/2)^2);%nm; height of the triangles
P.actin = unique(edges_s);
edge_type = zeros(size(edges_s,1),1);
P = check_free_T(T_s,edges_s,edge_type,P);
P.myosin_Tfree0 = P.myosin_Tfree;
% initialize mayosin edges
P.d_inter = 2;%auxiliar variable to check myosin intersections
P.max_r = 5*P.y0/2;%nm, maximum lenght of myosins rods
P.min_r = 3*P.y0/4;%nm, minimum lenght of myosins rods
P.myosin_T = []; %P.myosin_T: entries of the triangulation T_s corresponding to the triangles that have myosin rods
P.myosin_ini = 100;%number of initial myosin rodsX2

% create myosin rods
myosin = [];
while size(myosin,1) < P.myosin_ini
    [r_s,edges_s,edge_type,myosin,P] = create_myosin_3D(r_s,T_s,edges_s,edge_type,myosin,P);
end
% myosin: entries of the positin vector r_s corresponding to
% myosin rods nodes


% % create a plot to check if the mesh is correct
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
% 
% aux_e = find(edge_type == 2);
% for l=1:size(aux_e,1)
%     plot3([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
%        [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
%        [r_s(edges_s(aux_e(l),1),3),r_s(edges_s(aux_e(l),2),3)],'r')
% end
% axis('equal')
% plot3(r_s(P.actin,1),r_s(P.actin,2),r_s(P.actin,3),'ms')%plot actin short filaments 
% plot3(r_s(myosin,1),r_s(myosin,2),r_s(myosin,3),'g^')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% model parameters
P.k0 = 1;%pN/nm,spectrin spring constant
P.d00 = P.y0;%nm, spectrin resting length

P.theta0 = 0;%degrees, spontaneous curvature angle
P.k_bend =200*4.1;%nm*pN, average bending modules
P.k_A = 600*4.1/(2*P.d00^2);%area constraint
% parameters for surface area constraint 
T_aux_i1 = T_s(P.myosin_Tfree0,1);
T_aux_i2 = T_s(P.myosin_Tfree0,2);
T_aux_i3 = T_s(P.myosin_Tfree0,3);
s1 = r_s(T_aux_i1,:);
s2 = r_s(T_aux_i2,:);
s3 = r_s(T_aux_i3,:);
N = cross(s2-s1,s3-s1);
N_l = sqrt(dot(N,N,2));
P.A0 = (sum(N_l)./2);


P.gamma = 0.1071;%pN/nm,myosin2 cable constant

P.psi_a = 1/100;%1/s,myosin2 addition rate
P.psi_r = 1/160;%1/s,myosin2 removal rate

P.th = 5e-2;%pN, force threshold
P.d_restore = P.y0;%P.k0*P.d00/(P.k0+P.th);%



P.delta_t =0.002;% 0.002;%s, time step size
P.zeta = 1.25;%pN s/nm, drag


P.t_ini=0;%s, initial time
P.t_end = 360;%360;%s,end time
aux_t = P.t_ini:P.delta_t:P.t_end;


P.edges_s0= edges_s;
P.edge_type0 = edge_type;
P.T_s0 = T_s;

%  save variables
save_aux = 10;
r_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
edges_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
type_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
unbinding = zeros(P.t_end/(save_aux*P.delta_t)+1,1);
binding = zeros(P.t_end/(save_aux*P.delta_t)+1,1);
mr_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
T_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);


ll=1;
r_save{ll,1} = r_s;
edges_save{ll,1} = edges_s;
type_save{ll,1} = edge_type;
edges_remove{ll,1} = [];
unbinding(ll,1) = 0;
binding(ll,1) =  0;
mr_save{ll,1} = myosin;
T_save{ll,1} = T_s;

e_remove = [];
aux_unbinding = 0;
aux_binding = 0;


for k=2:length(aux_t)

    f = zeros(size(r_s));
    edges_s1 = [];
    aux_m = 1;
    edges_s2 = [];

   
    for l = 1:size(edges_s,1) 
%         calculate the distance of edges
        r_ij = r_s(edges_s(l,1),:) - r_s(edges_s(l,2),:);
        d = sqrt(dot(r_ij,r_ij,2));
%         depending of the type of edge compute the force and store the
%         value

%         for spectrin edges
        if edge_type(l) == 0
            if P.k0*(d-P.d00)/d > -P.th 
                f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.k0*r_ij*(d-P.d00)/d;
                f(edges_s(l,2),:) =  f(edges_s(l,2),:) + P.k0*r_ij*(d-P.d00)/d;
            else
                edges_s1 = [edges_s1;l];
            end
%         for myosin edges
        elseif edge_type(l) == 2
            aux_a = T_s(P.myosin_T(2*aux_m-1),:);
%             distribute the force among the three nodes of the pectrin
%             triangles that the myosin rods attaches to 
            f(aux_a,:) =  f(aux_a,:) - P.gamma*r_ij/3;
            aux_b = T_s(P.myosin_T(2*aux_m),:);
            f(aux_b,:) =  f(aux_b,:) + P.gamma.*r_ij/3;
            aux_m = aux_m + 1;
        end
    end
%     update the nodes position 
    f_d = force_def_area(r_s,P);%calculate force with area constratint 
    r_s(P.actin,:) = P.delta_t*(f(P.actin,:)+f_d(P.actin,:))./P.zeta + r_s(P.actin,:);
    
    
    %  for the spectrin edges that were detached, get the corresponding triangles  
    if ~isempty(edges_s1)
        for l = 1:length(edges_s1)
            aux_T = intersect([find(T_s(:,1) == edges_s(edges_s1(l),1));
                    find(T_s(:,2) == edges_s(edges_s1(l),1));
                    find(T_s(:,3) == edges_s(edges_s1(l),1))],...
                    [find(T_s(:,1) == edges_s(edges_s1(l),2));
                    find(T_s(:,2) == edges_s(edges_s1(l),2));
                    find(T_s(:,3) == edges_s(edges_s1(l),2))]);
            P.myosin_Tfree = setdiff(P.myosin_Tfree,aux_T);
        end
        
        for l = 1:length(edges_s1)
%             check if the myosin rod in the spectrin triangle that lost an
%             edge can attach to another spectrin triangle 
            [r_s,myosin,edges_s2,P] = check_myosin_3D(T_s,r_s,edges_s,edge_type,myosin,edges_s1(l),edges_s2,P);
        end
    end
    
    %     randomly add a myosin rod
    if rand <= P.psi_a*P.delta_t
       [r_s,edges_s,edge_type,myosin,P] = create_myosin_3D(r_s,T_s,edges_s,edge_type,myosin,P); 
    end
    
%     update position of the myosins

    for l = 1:length(myosin)
        r_s(myosin(l),:) = sum(r_s(T_s(P.myosin_T(l),:),:))./3;
    end
%     update unbound events
    aux_unbinding = aux_unbinding + size(edges_s1,1);
    e_remove = [e_remove;edges_s(edges_s1,:)];
    edges_s([edges_s1;edges_s2],:) = [];
    edge_type([edges_s1;edges_s2],:) = [];
    
    %     remove small myosin
    [edges_s,edge_type,myosin,P] = remove_small_myosin(r_s,edges_s,edge_type,myosin,P);
%     randomly remove a myosin rod
    if rand <= P.psi_r*P.delta_t
        [edges_s,edge_type,myosin,P] = remove_myosin(edges_s,edge_type,myosin,P);
    end
%     check if unbound spectrin edges can be re-bind 
    if ~isempty(e_remove)
        [e_remove,edges_s,edge_type,aux_binding,T_s,P] = reatach_spectrin(e_remove,edges_s,edge_type,aux_binding,r_s,T_s,P);
    end
    
%     save data 
    if mod(k,save_aux) == 1
        
        ll=ll+1;
        r_save{ll,1} = r_s;
        edges_save{ll,1} = edges_s;
        type_save{ll,1} = edge_type;
        edges_remove{ll,1} = e_remove;
        mr_save{ll,1} = myosin;
        T_save{ll,1} = T_s;
        unbinding(ll,1) = aux_unbinding;
        binding(ll,1) =  aux_binding;
        aux_unbinding = 0;
        aux_binding = 0;
        
    end
end
