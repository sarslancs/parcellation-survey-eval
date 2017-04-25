function [ silhouettes ] = silhouette_coef( parcels, D, neigh )
%SILHOUETTE_COEF Silhouette coefficient of a parcellation.
%   For each vertex in a cortical surface, SILHOUETTE_COEF compares the 
%   within-parcel dissimilarity defined as the average distance to all 
%   vertices in the same parcel, to the inter-parcel dissimilarity 
%   obtained from those assigned to other parcels.
%
%   Let Pi be the parcel to which vertex i is assigned and ai and bi be the
%   average distance from i to the vertices in Pi and to the vertices in 
%   other parcels adjacent to Pi. Silhouette coefficient (SI) for i is then 
%   computed as:
%
%   SI(i) = (bi - ai)/max(ai,bi) 
%
%   This guarantees SI values within [-1, +1], as long as a distance 
%   measure is used. 
%
%   INPUT
%   =====
%   parcels: A parcellation.
%   D: A distance matrix (such as Pearson's distance).
%   neigh: An adjacency matrix that provides the "neighbourhood" 
%          information for each parcel and is used to identify the vertices
%          in adjacent parcels. Two parcels A and B are considered as 
%          neighbours, or adjacent, if vertices i ∈ A and j ∈ B are 
%          directly connected by an edge in the cortical mesh. Such
%          matrices are provided for each parcellation.
%
%   OUTPUT
%   ======
%   silhouettes: Silhouette coefficients for all vertices
%
%   USAGE
%   =====
%   [ SILS ] = SILHOUETTE_COEF( PARCELS, D, NEIGH ) returns an N-by-1 
%   vector, in which the ith element indicates the Silhouette value of the 
%   ith vertex and N is the number of vertices. PARCELS can be a 
%   parcellation of any resolution. D must be an N-by-N matrix, where each
%   vertex pair (x,y) equals to the distance between x and y. NEIGH is a
%   K-by-K adjaceny matrix, where K denotes the parcellation resolution.
%
%   REFERENCE
%   This code is part of the evaluation pipelines described in the brain
%   parcellation survey, "Human Brain Mapping: A Systematic Comparison of
%   Parcellation Methods for the Human Cerebral Cortex", NeuroImage, 2017
%   doi.org/10.1016/j.neuroimage.2017.04.014 
%
%   For parcellations and more visit the Brain Parcellation Survey page at 
%   https://biomedia.doc.ic.ac.uk/brain-parcellation-survey/ 
%
%   Author: Salim Arslan, April 2017 (twitter @salimarslan)

n = length(D);
num = max(parcels);
silhouettes = zeros(n,1);

for i = 1 : num  
    in_members = parcels == i;
    nk = sum(in_members); 
    
    if nk < 2 
        continue; % Singleton parcel detected. Can happen with N-Cuts.
    end
        
    dists = D(in_members,in_members);
    dists(logical(eye(length(dists)))) = 0;
    dists_in = sum(dists,2)/(length(dists)-1); 
  
    ids = find(neigh(:,i));
    dists_out = find_dist_to_neighs( ids, parcels, D, in_members );
    
    silhouettes(in_members) = (dists_out - dists_in) ./ ...
                              max([dists_in dists_out],[],2); 
   
end

assert(sum(isnan(silhouettes))==0);

end


% Compute distance to/from vertices in adjacent parcels
function [ votes_out] = find_dist_to_neighs(ids, parcels, D, in_members )

out_members = false(size(parcels));
for j = 1 : length(ids)
    out_members = out_members | parcels == ids(j);
end
nk = sum(out_members); 
corrs = D(in_members,out_members);
votes_out = sum(corrs,2)/nk;

end





