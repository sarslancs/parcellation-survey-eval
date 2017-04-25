function [ homogeneities ] = homogeneity( parcels, Z )
%HOMOGENEITY Parcel homogeneity.
%   Homogeneity of a parcel is measured by calculating the average
%   similarity between every pair of vertices assigned to it. A global 
%   homogeneity value for the entire parcellation can be later obtained by 
%   averaging the homogeneity values across all parcels.
%
%   INPUT
%   =====
%   parcels: A parcellation.
%   Z: A correlation matrix (ideally Fisher's r-to-z transformed).
%
%   OUTPUT
%   ======
%   homogeneities: Homogeneity values.
%
%   USAGE
%   =====
%   HOMS = HOMOGENEITY( PARCELS, Z ) returns a K-by-1 vector, in which the 
%   kth element indicates the homogeneity value of the kth parcel and K is
%   the number of labels. PARCELS is an N-by-1 parcellation vector, where 
%   N denotes the number of vertices. Z must be an N-by-N matrix, in which
%   each vertex pair (x,y) equals to the correlation of x and y. 
%
%   REFERENCE
%   =========
%   This code is part of the evaluation pipelines described in the brain
%   parcellation survey, "Human Brain Mapping: A Systematic Comparison of
%   Parcellation Methods for the Human Cerebral Cortex", NeuroImage, 2017
%   doi.org/10.1016/j.neuroimage.2017.04.014 
%
%   For parcellations and more visit the Brain Parcellation Survey page at 
%   https://biomedia.doc.ic.ac.uk/brain-parcellation-survey/ 
%
%   Author: Salim Arslan, April 2017 (name.surname@imperial.ac.uk)

K = max(parcels);
homogeneities = zeros(K,1);
for i = 1 : K  
    in_members = parcels == i;
    nk = sum(in_members); 
    
    if nk < 2 % In case there are parcels with only 1 element (may happen with NCUTS)
        ak = 1;
    else
        corrs = Z(in_members,in_members)';
        corrs(logical(eye(length(corrs)))) = 0;
        means_in = sum(corrs,2)/(nk-1);
        ak = mean(means_in);
    end
    homogeneities(i) = ak;
end