function [r_s,edges_s,edge_type,myosin,P] = create_myosin_3D(r_s,T_s,edges_s,edge_type,myosin,P)
%     select one free triangle and check the distance to nearby triangles
    aux_r1 = randi(length(P.myosin_Tfree));
    center_Tfree = (r_s(T_s(P.myosin_Tfree,1),:) + r_s(T_s(P.myosin_Tfree,2),:) +...
        r_s(T_s(P.myosin_Tfree,3),:))./3;
    aux_d = center_Tfree - center_Tfree(aux_r1,:);
    aux_d = sqrt(dot(aux_d,aux_d,2));
    aux_new = find(aux_d < P.max_r & aux_d > P.min_r);
    
    if ~isempty(aux_new)
        
        aux_c = 1:length(aux_new);
        aux_r = randi(length(aux_c));
        aux = find(edge_type == 2);

        
        if ~isempty(aux) 
% check if there are myosin rods inside the sphere generated by the new
% rods
            a = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
            b = center_Tfree(aux_r1,:);
            
            p = r_s(edges_s(aux,1),:);
            q = r_s(edges_s(aux,2),:);

            r = sqrt(dot(a-b,a-b,2))./2;
            c = (a+b)./2;
            
            cond1 = find(sqrt(dot(p-c,p-c,2))< P.d_inter*r ,1);
            cond2 = find(sqrt(dot(q-c,q-c,2))< P.d_inter*r ,1);
            
            while (~isempty(cond1) || ~isempty(cond2))  && length(aux_c)>1 
              
                aux_c(aux_r) = [];
                aux_r = randi(length(aux_c));
                
                a = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
                r = sqrt(dot(a-b,a-b,2))./2;
                c = (a+b)./2;
            
                cond1 = find(sqrt(dot(p-c,p-c,2))< P.d_inter*r ,1);
                cond2 = find(sqrt(dot(q-c,q-c,2))< P.d_inter*r ,1);
                
                
            end

            if isempty(cond1) && isempty(cond2) 
                P.myosin_T =  [P.myosin_T;P.myosin_Tfree(aux_r1);P.myosin_Tfree(aux_new(aux_c(aux_r)))];
                P.myosin_Tfree([aux_new(aux_c(aux_r));aux_r1]) = [];
                edges_s = [edges_s;size(r_s,1)+(1:2)];
                edge_type = [edge_type;2];
                myosin = [myosin;size(r_s,1)+(1:2)'];
                r_s = [r_s;sum(r_s(T_s(P.myosin_T(end-1),:),:))./3;sum(r_s(T_s(P.myosin_T(end),:),:))./3];
            end
        else
            aux_r = randi(length(aux_new));
            P.myosin_T =  [P.myosin_T;P.myosin_Tfree(aux_r1);P.myosin_Tfree(aux_new(aux_r))];
            P.myosin_Tfree([aux_new(aux_r);aux_r1]) = [];
            edges_s = [edges_s;size(r_s,1)+(1:2)];
            edge_type = [edge_type;2];
            myosin = [myosin;size(r_s,1)+(1:2)'];
            r_s = [r_s;sum(r_s(T_s(P.myosin_T(end-1),:),:))./3;sum(r_s(T_s(P.myosin_T(end),:),:))./3];
        end
    end
end