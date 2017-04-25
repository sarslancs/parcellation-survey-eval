function [ relabeled, K ] = relabel( parcels )
%RELABEL Relabel a parcellation. Just a utility function.

ids = nonzeros(unique((parcels)));
K = length(ids);
if max(parcels) == K
    relabeled = parcels;
else
    relabeled = zeros(size(parcels));
    id = 1;
    for i = 1 : K
        relabeled(parcels == ids(i)) = id;
        id = id + 1;
    end
end





