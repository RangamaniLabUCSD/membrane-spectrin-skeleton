function [edges_s,edge_type,myosin,P] = remove_small_myosin(r_s,edges_s,edge_type,myosin,P)
% check the length of the myosin rods
    aux = find(edge_type == 2);
    d = r_s(edges_s(aux,1),:) - r_s(edges_s(aux,2),:);
    d = sqrt(dot(d,d,2));
    aux_r = find(d < P.min_r);
    %     remove those myosin rods that are shorter than the minimum length
    if ~isempty(aux_r)
        aux_rem = [];
        for l = 1:length(aux_r)
            aux_r1 = find(myosin == edges_s(aux(aux_r(l)),1));
            if ~isempty(aux_r1)
                aux_r2 = find(myosin == edges_s(aux(aux_r(l)),2));
            else
                aux_r1 = find(myosin == edges_s(aux(aux_r(l)),2));
                aux_r2 = find(myosin == edges_s(aux(aux_r(l)),1));
            end
            aux_rem = [aux_rem;aux_r1;aux_r2];
        end
        myosin(aux_rem) = [];
        P.myosin_Tfree = [P.myosin_Tfree;P.myosin_T(aux_rem)];
        P.myosin_T(aux_rem) = [];
        edges_s(aux(aux_r),:) = [];
        edge_type(aux(aux_r)) = [];
    end
end