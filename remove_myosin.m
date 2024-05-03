function [edges_s,edge_type,myosin,P] = remove_myosin(edges_s,edge_type,myosin,P)
% randomly remove one myosin rod 
    aux = find(edge_type == 2);
    if ~isempty(aux)
        aux_e = randi(length(aux));
        aux_r1 = find(myosin == edges_s(aux(aux_e),1));
        if ~isempty(aux_r1)
            aux_r2 = find(myosin == edges_s(aux(aux_e),2));
        else
            aux_r1 = find(myosin == edges_s(aux(aux_e),2));
            aux_r2 = find(myosin == edges_s(aux(aux_e),1));
        end
        myosin([aux_r1;aux_r2]) = [];
        P.myosin_Tfree = [P.myosin_Tfree;P.myosin_T([aux_r1;aux_r2])];
        P.myosin_T([aux_r1;aux_r2]) = [];
        edges_s(aux(aux_e),:) = [];
        edge_type(aux(aux_e)) = [];
    end
end