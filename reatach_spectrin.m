function [e_remove,edges_s,edge_type,aux_binding,T_s,P] = reatach_spectrin(e_remove,edges_s,edge_type,aux_binding,r_s,T_s,P)
%     from the unbound spectrin check which can be re-bind 
    r_ij = r_s(e_remove(:,1),:) - r_s(e_remove(:,2),:);
    d = sqrt(dot(r_ij,r_ij,2));
    edges_change = e_remove(d>P.d_restore,:);
    e_remove(d>P.d_restore,:) = [];
    edges_s = [edges_s;edges_change];
    edge_type = [edge_type;zeros(size(edges_change,1),1)];
    aux_binding = aux_binding + size(edges_change,1);
    
%     add free triangle 
    edges0 = edges_s(edge_type == 0,:);
    e3 = [];
    e3b = [];
    e4 = [];
    e4b = [];
    
    
  
    for l = 1:size(edges_change,1)
        e1 = find(edges0(:,1) == edges_change(l,1));
        if ~isempty(e1)
            for ll=1:length(e1)
                e3 = [e3;find(edges0(:,2) == e1(ll))];
            end
        end
        e1b = find(edges0(:,2) == edges_change(l,1));
        if ~isempty(e1b)
            for ll=1:length(e1b)
                e3b = [e3b;find(edges0(:,1) == e1b(ll))];
            end
        end
        e2 = find(edges0(:,1) == edges_change(l,2));
        if ~isempty(e2)
            for ll=1:length(e2)
                e4 = [e4;find(edges0(:,2) == e2(ll))];
            end
        end
        e2b = find(edges0(:,2) == edges_change(l,2));
        if ~isempty(e2b)
            for ll=1:length(e2b)
                e4b = [e4b;find(edges0(:,1) == e2b(ll))];
            end
        end
        
        aux = intersect(edges0(e3,2),edges0(e4,2));
        if ~isempty(aux)
            for ll=1:length(aux)
                T_s = [T_s;edges_change(l,1) edges_change(l,2) aux(ll)];
            end
        end
        aux = intersect(edges0(e3,2),edges0(e4b,1));
        if ~isempty(aux)
            for ll=1:length(aux)
                T_s = [T_s;edges_change(l,1) edges_change(l,2) aux(ll)];
            end
        end
        
        aux = intersect(edges0(e3b,1),edges0(e4,2));
        if ~isempty(aux)
            for ll=1:length(aux)
                T_s = [T_s;edges_change(l,1) edges_change(l,2) aux(ll)];
            end
        end
        aux = intersect(edges0(e3b,1),edges0(e4b,1));
        if ~isempty(aux)
            for ll=1:length(aux)
                T_s = [T_s;edges_change(l,1) edges_change(l,2) aux(ll)];
            end
        end
    end


end