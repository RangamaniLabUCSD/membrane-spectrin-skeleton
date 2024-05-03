function f_d = force_def(S,P)
    edges_s = P.edges_s0;
    edge_type = P.edge_type0;
    T_s = P.T_s0;
    
    f_d = zeros(size(S));%save force generated by bending
    A_av = [];%save triangle area
    
    for k = 1:size(edges_s,1)
        if edge_type(k) == 0
            aux1 = edges_s(k,1);
            aux2 = edges_s(k,2);
%             define the position of the edges in the triangle 
            if ismember(aux1,P.actin) &&  ismember(aux2,P.actin)
                si = [];
                sip1 = [];
                sip2 = [];
                sj = [];
                sjp1 = [];
                sjp2 = [];


                aux_Ti = find(T_s(:,1)==aux1);
                aux_Ti = intersect(P.myosin_Tfree0,aux_Ti);
                if ~isempty(aux_Ti)
                    aux_Tip1 = find(T_s(aux_Ti,2)==aux2);
                    if ~isempty(aux_Tip1)
                        si = T_s(aux_Ti(aux_Tip1),1);
                        sip1 = T_s(aux_Ti(aux_Tip1),2);
                        sip2 = T_s(aux_Ti(aux_Tip1),3);
                        sj = sip1;
                        sjp1 = si;

                        aux_j = find(T_s(:,1)==aux2);
                        aux_j = intersect(P.myosin_Tfree0,aux_j);
                        if ~isempty(aux_j)
                            aux_jp = find(T_s(aux_j,2)==aux1);
                            if ~isempty(aux_jp)
                                sjp2 = T_s(aux_j(aux_jp),3);
                            end
                        else
                            aux_j = find(T_s(:,2)==aux2);
                            aux_j = intersect(P.myosin_Tfree0,aux_j);
                            if ~isempty(aux_j)
                                aux_jp = find(T_s(aux_j,3)==aux1);
                                if ~isempty(aux_jp)
                                    sjp2 = T_s(aux_j(aux_jp),1);
                                end
                            else
                                aux_j = find(T_s(:,3)==aux2);
                                aux_j = intersect(P.myosin_Tfree0,aux_j);
                                if ~isempty(aux_j)
                                    aux_jp = find(T_s(aux_j,1)==aux1);
                                    if ~isempty(aux_jp)
                                        sjp2 = T_s(aux_j(aux_jp),2);
                                    end
                                end
                            end
                        end
                    else
                        aux_Ti = find(T_s(:,2)==aux1);
                        aux_Ti = intersect(P.myosin_Tfree0,aux_Ti);
                        if ~isempty(aux_Ti)
                            aux_Tip1 = find(T_s(aux_Ti,3)==aux2);
                            if ~isempty(aux_Tip1)
                                si = T_s(aux_Ti(aux_Tip1),2);
                                sip1 = T_s(aux_Ti(aux_Tip1),3);
                                sip2 = T_s(aux_Ti(aux_Tip1),1);
                                sj = sip1;
                                sjp1 = si;

                                aux_j = find(T_s(:,1)==aux2);
                                aux_j = intersect(P.myosin_Tfree0,aux_j);
                                if ~isempty(aux_j)
                                    aux_jp = find(T_s(aux_j,2)==aux1);
                                    if ~isempty(aux_jp)
                                        sjp2 = T_s(aux_j(aux_jp),3);
                                    end
                                else
                                    aux_j = find(T_s(:,2)==aux2);
                                    aux_j = intersect(P.myosin_Tfree0,aux_j);
                                    if ~isempty(aux_j)
                                        aux_jp = find(T_s(aux_j,3)==aux1);
                                        if ~isempty(aux_jp)
                                            sjp2 = T_s(aux_j(aux_jp),1);
                                        end
                                    else
                                        aux_j = find(T_s(:,3)==aux2);
                                        aux_j = intersect(P.myosin_Tfree0,aux_j);
                                        if ~isempty(aux_j)
                                            aux_jp = find(T_s(aux_j,1)==aux1);
                                            if ~isempty(aux_jp)
                                                sjp2 = T_s(aux_j(aux_jp),2);
                                            end
                                        end
                                    end
                                end
                            else
                                aux_Ti = find(T_s(:,3)==aux1);
                                aux_Ti = intersect(P.myosin_Tfree0,aux_Ti);
                                if ~isempty(aux_Ti)
                                    aux_Tip1 = find(T_s(aux_Ti,1)==aux2);
                                    if ~isempty(aux_Tip1)
                                        si = T_s(aux_Ti(aux_Tip1),3);
                                        sip1 = T_s(aux_Ti(aux_Tip1),1);
                                        sip2 = T_s(aux_Ti(aux_Tip1),2);
                                        sj = sip1;
                                        sjp1 = si;

                                        aux_j = find(T_s(:,1)==aux2);
                                        aux_j = intersect(P.myosin_Tfree0,aux_j);
                                        if ~isempty(aux_j)
                                            aux_jp = find(T_s(aux_j,2)==aux1);
                                            if ~isempty(aux_jp)
                                                sjp2 = T_s(aux_j(aux_jp),3);
                                            end
                                        else
                                            aux_j = find(T_s(:,2)==aux2);
                                            aux_j = intersect(P.myosin_Tfree0,aux_j);
                                            if ~isempty(aux_j)
                                                aux_jp = find(T_s(aux_j,3)==aux1);
                                                if ~isempty(aux_jp)
                                                    sjp2 = T_s(aux_j(aux_jp),1);
                                                end
                                            else
                                                aux_j = find(T_s(:,3)==aux2);
                                                aux_j = intersect(P.myosin_Tfree0,aux_j);
                                                if ~isempty(aux_j)
                                                    aux_jp = find(T_s(aux_j,1)==aux1);
                                                    if ~isempty(aux_jp)
                                                        sjp2 = T_s(aux_j(aux_jp),2);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                if isempty(si)
                    aux1 = edges_s(k,2);
                    aux2 = edges_s(k,1);
                    aux_Tj = find(T_s(:,1)==aux1);
                    aux_Tj = intersect(P.myosin_Tfree0,aux_Tj);
                    if ~isempty(aux_Tj)
                        aux_Tjp1 = find(T_s(aux_Tj,2)==aux2);
                        if ~isempty(aux_Tjp1)
                            sj = T_s(aux_Tj(aux_Tjp1),1);
                            sjp1 = T_s(aux_Tj(aux_Tjp1),2);
                            sjp2 = T_s(aux_Tj(aux_Tjp1),3);
                            si = sjp1;
                            sip1 = sj;

                            aux_i = find(T_s(:,1)==aux2);
                            aux_i = intersect(P.myosin_Tfree0,aux_i);
                            if ~isempty(aux_i)
                                aux_ip = find(T_s(aux_i,2)==aux1);
                                if ~isempty(aux_ip)
                                    sip2 = T_s(aux_i(aux_ip),3);
                                end
                            else
                                aux_i = find(T_s(:,2)==aux2);
                                aux_i = intersect(P.myosin_Tfree0,aux_i);
                                if ~isempty(aux_i)
                                    aux_ip = find(T_s(aux_i,3)==aux1);
                                    if ~isempty(aux_ip)
                                        sip2 = T_s(aux_i(aux_ip),1);
                                    end
                                else
                                    aux_i = find(T_s(:,3)==aux2);
                                    aux_i = intersect(P.myosin_Tfree0,aux_i);
                                    if ~isempty(aux_i)
                                        aux_ip = find(T_s(aux_i,1)==aux1);
                                        if ~isempty(aux_ip)
                                            sip2 = T_s(aux_i(aux_ip),2);
                                        end
                                    end
                                end
                            end
                        else
                            aux_Tj = find(T_s(:,2)==aux1);
                            aux_Tj = intersect(P.myosin_Tfree0,aux_Tj);
                            if ~isempty(aux_Tj)
                                aux_Tjp1 = find(T_s(aux_Tj,3)==aux2);
                                if ~isempty(aux_Tjp1)
                                    sj = T_s(aux_Tj(aux_Tjp1),2);
                                    sjp1 = T_s(aux_Tj(aux_Tjp1),3);
                                    sjp2 = T_s(aux_Tj(aux_Tjp1),1);
                                    si = sjp1;
                                    sip1 = sj;

                                    aux_i = find(T_s(:,1)==aux2);
                                    aux_i = intersect(P.myosin_Tfree0,aux_i);
                                    if ~isempty(aux_i)
                                        aux_ip = find(T_s(aux_i,2)==aux1);
                                        if ~isempty(aux_ip)
                                            sip2 = T_s(aux_i(aux_ip),3);
                                        end
                                    else
                                        aux_i = find(T_s(:,2)==aux2);
                                        aux_i = intersect(P.myosin_Tfree0,aux_i);
                                        if ~isempty(aux_i)
                                            aux_ip = find(T_s(aux_i,3)==aux1);
                                            if ~isempty(aux_ip)
                                                sip2 = T_s(aux_i(aux_ip),1);
                                            end
                                        else
                                            aux_i = find(T_s(:,3)==aux2);
                                            aux_i = intersect(P.myosin_Tfree0,aux_i);
                                            if ~isempty(aux_i)
                                                aux_ip = find(T_s(aux_i,1)==aux1);
                                                if ~isempty(aux_ip)
                                                    sip2 = T_s(aux_i(aux_ip),2);
                                                end
                                            end
                                        end
                                    end
                                else
                                    aux_Tj = find(T_s(:,3)==aux1);
                                    aux_Tj = intersect(P.myosin_Tfree0,aux_Tj);
                                    if ~isempty(aux_Tj)
                                        aux_Tjp1 = find(T_s(aux_Tj,1)==aux2);
                                        if ~isempty(aux_Tjp1)
                                            sj = T_s(aux_Tj(aux_Tjp1),3);
                                            sjp1 = T_s(aux_Tj(aux_Tjp1),2);
                                            sjp2 = T_s(aux_Tj(aux_Tjp1),1);
                                            si = sjp1;
                                            sip1 = sj;

                                            aux_i = find(T_s(:,1)==aux2);
                                            aux_i = intersect(P.myosin_Tfree0,aux_i);
                                            if ~isempty(aux_i)
                                                aux_ip = find(T_s(aux_i,2)==aux1);
                                                if ~isempty(aux_ip)
                                                    sip2 = T_s(aux_i(aux_ip),3);
                                                end
                                            else
                                                aux_i = find(T_s(:,2)==aux2);
                                                aux_i = intersect(P.myosin_Tfree0,aux_i);
                                                if ~isempty(aux_i)
                                                    aux_ip = find(T_s(aux_i,3)==aux1);
                                                    if ~isempty(aux_ip)
                                                        sip2 = T_s(aux_i(aux_ip),1);
                                                    end
                                                else
                                                    aux_i = find(T_s(:,3)==aux2);
                                                    aux_i = intersect(P.myosin_Tfree0,aux_i);
                                                    if ~isempty(aux_i)
                                                        aux_ip = find(T_s(aux_i,1)==aux1);
                                                        if ~isempty(aux_ip)
                                                            sip2 = T_s(aux_i(aux_ip),2);
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end

                end
%                 if it is a complete triangle calculate the force
%                 generated by the bending energy 
                if  ~isempty(sjp2) && ~isempty(sip2)
                    if ismember(sjp2,P.actin) &&  ismember(sip2,P.actin)
                        Ni = cross(S(sip1,:)-S(si,:),S(sip2,:)-S(si,:));
                        length_Ni = sqrt(dot(Ni,Ni,2));
                        ni = Ni./length_Ni;
                        Nj = cross(S(sjp1,:)-S(sj,:),S(sjp2,:)-S(sj,:));
                        length_Nj = sqrt(dot(Nj,Nj,2));
                        nj = Nj./length_Nj;
                        Aij = sqrt(dot(Ni,Ni,2))*sqrt(dot(Nj,Nj,2))./4;
                        A_av =[A_av;Aij];

                        i_c = (S(si,:)+S(sip1,:)+S(sip2,:))./3;
                        j_c = (S(sj,:)+S(sjp1,:)+S(sjp2,:))./3;
                        aux_sign = dot(ni-nj,i_c-j_c,2);
                        if aux_sign >=0
                            aux_sign = 1;
                        else
                            aux_sign = -1;
                        end
                        aux_n = dot(ni,nj,2);

                        aux_cross = cross(ni,nj);
                        length_cross = dot(aux_cross,aux_cross,2);
                        if length_cross == 0
                            sin_si = [0 0 0];
                            sin_sip1 = [0 0 0];
                            sin_sip2 = [0 0 0];
                            sin_sjp2 = [0 0 0];
                        else
                            sin_si_i = -cross(ni,S(sip2,:)-S(sip1,:)).*aux_cross./length_Ni;
                            sin_si_j = -aux_cross.*cross(nj,S(sip1,:)-S(sjp2,:))./length_Nj;
                            sin_si = aux_sign*dot(aux_cross,sin_si_i+sin_si_j,2)./length_cross;
        %     
                            sin_sip1_i = -cross(S(sip2,:)-S(si,:),ni).*aux_cross./length_Ni;
                            sin_sip1_j = -aux_cross.*cross(S(si,:)-S(sjp2,:),nj)./length_Nj;
                            sin_sip1 = aux_sign*dot(aux_cross,sin_sip1_i+sin_sip1_j,2)./length_cross;

                            sin_sip2 =  aux_sign*dot(aux_cross,-aux_cross.*cross(ni,S(sip1,:)-S(si,:))./length_Ni,2)./length_cross;
                            sin_sjp2 = aux_sign*dot(aux_cross,-aux_cross.*cross(S(sip1,:)-S(si,:),nj)./length_Nj,2)./length_cross;
                        end


%                         update the forces for each node
                        f_d(si,:) = f_d(si,:) + Aij.*(cos(P.theta0)*(cross(nj-aux_n.*ni,S(sip2,:)-S(sip1,:))./length_Ni ...
                            + cross(ni-aux_n.*nj,S(sj,:)-S(sjp2,:))./length_Nj)...
                            +sin(P.theta0).*sin_si);

                        f_d(sip1,:) = f_d(sip1,:) + Aij.*(cos(P.theta0)*(cross(nj-aux_n.*ni,S(si,:)-S(sip2,:))./length_Ni ...
                            + cross(ni-aux_n.*nj,S(sjp2,:)-S(sjp1,:))./length_Nj)...
                            +sin(P.theta0).*sin_sip1);

                        f_d(sip2,:) = f_d(sip2,:) + Aij.*(cos(P.theta0)*(cross(nj-aux_n.*ni,S(sip1,:)-S(si,:))./length_Ni)...
                            +sin(P.theta0).*sin_sip2);

                        f_d(sjp2,:) = f_d(sjp2,:) + Aij.*(cos(P.theta0)*(cross(ni-aux_n.*nj,S(sjp1,:)-S(sj,:))./length_Nj)...
                            +sin(P.theta0).*sin_sjp2);
                    end
                end
            end
        end
    end
    f_d = P.k_bend*f_d./mean(A_av);
   
end
