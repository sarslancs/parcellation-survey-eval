function [ ari ] = adj_rand_index(U,V)
%ADJ_RAND_INDEX Adjusted rand index (ARI).
%   Computes the adjusted Rand index to assess the similarity of 
%   parcellations in a scan-to-scan or group-to-group setting.
%   An ARI of 0 indicates perfectly random clustering, while equivalent
%   parcellations return the maximum score of 1.
%
%   INPUT
%   =====
%   U: First parcellation, labeled from 1 to K1 
%   V: Second parcellation, labeled from 1 to K2 (K1 is not necessarily 
%      equal to K2)
%
%   OUTPUT
%   ======
%   adjrand: the adjusted Rand index 
%
%   USAGE
%   =====
%   ARI = ADJ_RAND_INDEX( U, V ) returns the ARI, which indicates the 
%   similarity between parcellations U and V. U and V can be (1) 
%   parcellations of the same subject acquired from different scans or 
%   (2) groupwise parcellations of a population. Input parcellations can
%   be of different resolution, but must be the same size.
%
%   REFERENCE
%   =========
%   This code is part of the evaluation pipelines described in the brain
%   parcellation survey, "Human Brain Mapping: A Systematic Comparison of
%   Parcellation Methods for the Human Cerebral Cortex", NeuroImage, 2017
%   doi.org/10.1016/j.neuroimage.2017.04.014 
%
%   For the parcellation data and reference manual visit the survey page: 
%   https://biomedia.doc.ic.ac.uk/brain-parcellation-survey/ 
%
%   Modified by Salim Arslan, April 2017 (name.surname@imperial.ac.uk)
%
%   Author: Tijl De Bie, February 2003 


n=length(U);
ku=max(U);
kv=max(V);
m=zeros(ku,kv);
for i=1:n
    m(U(i),V(i))=m(U(i),V(i))+1;
end
mu=sum(m,2);
mv=sum(m,1);

a=0;
for i=1:ku
    for j=1:kv
        if m(i,j)>1
            a=a+nchoosek(m(i,j),2);
        end
    end
end

b1=0;
b2=0;
for i=1:ku
    if mu(i)>1
        b1=b1+nchoosek(mu(i),2);
    end
end
for i=1:kv
    if mv(i)>1
        b2=b2+nchoosek(mv(i),2);
    end
end

c=nchoosek(n,2);

ari=(a-b1*b2/c)/(0.5*(b1+b2)-b1*b2/c);