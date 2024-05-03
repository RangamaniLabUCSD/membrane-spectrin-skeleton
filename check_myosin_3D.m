function [r_s,myosin,edges_s2,P] = check_myosin_3D(T_s,r_s,edges_s,edge_type,myosin,edges_s1,edges_s2,P)
%     check if detached myosin edges can attach to a nearby triangle
    aux_T1 = [find(T_s(P.myosin_T,1) == edges_s(edges_s1,1));
        find(T_s(P.myosin_T,2) == edges_s(edges_s1,1));
        find(T_s(P.myosin_T,3) == edges_s(edges_s1,1))];
    aux_T2 = [find(T_s(P.myosin_T,1) == edges_s(edges_s1,2));
        find(T_s(P.myosin_T,2) == edges_s(edges_s1,2));
        find(T_s(P.myosin_T,3) == edges_s(edges_s1,2))];
    aux_T3 = intersect(aux_T1,aux_T2);
    if ~isempty(aux_T3)
        aux_TT = [];
        aux_TT2 = [];
%         check distance to other free triangles
        for lll = 1:length(aux_T3)
            center_Tfree = (r_s(T_s(P.myosin_Tfree,1),:) + r_s(T_s(P.myosin_Tfree,2),:) +...
                r_s(T_s(P.myosin_Tfree,3),:))./3;
            aux = find(edge_type == 2);
            aux2 = [find(edges_s(aux,1) == myosin(aux_T3(lll)));
                    find(edges_s(aux,2) == myosin(aux_T3(lll)))];
            aux3 = [find(myosin == edges_s(aux(aux2),1));
                find(myosin == edges_s(aux(aux2),2))];
            aux4 = setdiff(aux3,aux_T3(lll));

            aux_d = center_Tfree - r_s(myosin(aux4),:);
            aux_d = sqrt(dot(aux_d,aux_d,2));
            aux_new = find(aux_d < P.max_r & aux_d > P.min_r);


            if ~isempty(aux_new)

                aux_c = 1:length(aux_new);
                aux_r = randi(length(aux_c));
                
                a = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
                b = r_s(myosin(aux4),:);

                p = r_s(edges_s(aux,1),:);
                q = r_s(edges_s(aux,2),:);

                r = sqrt(dot(a-b,a-b,2))./2;
                c = (a+b)./2;

                cond1 = find(sqrt(dot(p-c,p-c,2))< P.d_inter*r ,1);
                cond2 = find(sqrt(dot(q-c,q-c,2))< P.d_inter*r ,1);
                
                while  (~isempty(cond1) || ~isempty(cond2))  && length(aux_c)>1 
                    
                    aux_c(aux_r) = [];
                    aux_r = randi(length(aux_c));

                    a = sum(r_s(T_s(P.myosin_Tfree(aux_new(aux_c(aux_r))),:),:))./3;
                    r = sqrt(dot(a-b,a-b,2))./2;
                    c = (a+b)./2;

                    cond1 = find(sqrt(dot(p-c,p-c,2))< P.d_inter*r ,1);
                    cond2 = find(sqrt(dot(q-c,q-c,2))< P.d_inter*r ,1);   

                end
                if isempty(cond1) && isempty(cond2) 
                    P.myosin_T(aux_T3(lll)) = P.myosin_Tfree(aux_new(aux_c(aux_r)));
                    r_s(myosin(aux_T3(lll)),:) = a;
                    P.myosin_Tfree(aux_new(aux_c(aux_r))) = [];
                else
                    edges_s2 = [edges_s2;aux(aux2)];
                    aux_TT = [aux_TT;aux3];
                    aux_TT2 = [aux_TT2;aux4];
                end
            else
%                 disp(aux(aux2))
                edges_s2 = [edges_s2;aux(aux2)];
                aux_TT = [aux_TT;aux3];
                aux_TT2 = [aux_TT2;aux4];
            end
        end
        P.myosin_Tfree = [P.myosin_Tfree;P.myosin_T(aux_TT2)];
        myosin(aux_TT) = [];
        P.myosin_T(aux_TT) = [];
    end
end