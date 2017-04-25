function [ dice, Umatched, Vmatched ] = dice_coef( U, V )
%DICE_COEF Dice coefficient of two parcellations.
%   Computes the Dice score to assess the similarity of parcellations in a
%   scan-to-scan or group-to-group setting. An iterative procedure is 
%   performed to indetify matching parcels before computing pairwise 
%   Dice overlaps. A global Dice score is obtained via averaging.
%
%   A Dice score of 0 indicates perfectly random clustering, while 
%   equivalent parcellations return the maximum score of 1.
%
%   INPUT
%   =====
%   U: First parcellation, labeled from 1 to K1 
%   V: Second parcellation, labeled from 1 to K2 (K1 is not necessarily 
%      equal to K2)
%
%   OUTPUT
%   ======
%   dice: Average Dice coefficient
%   Umatched: Parcellation U, with new labels after being matched to V
%   Vmatched: Parcellation V, with new labels after being matched to U
%   
%   USAGE 
%   =====
%   DICE = DICE_COEF( U, V ) returns a DICE value, which indicates the 
%   similarity between parcellations U and V. U and V can be (1) 
%   parcellations of the same subject acquired from different scans or 
%   (2) groupwise parcellations of the same population. Input 
%   parcellations can be of different resolution, but must be the same 
%   size.
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
%   Author: Salim Arslan, April 2017 (name.surname@imperial.ac.uk)


% Relabel U and V in case of any missing labels 
[ parcels1, K1 ] = relabel( U );
[ parcels2, K2 ] = relabel( V );

Umatched = zeros(size(parcels1));
Vmatched = zeros(size(parcels2));

check1 = zeros(K1, 1);
check2 = zeros(K2, 1);

overlaps = zeros(max(K1,K2),1);
pairs = zeros(max(K1,K2),2);
dices = zeros(max(K1,K2),1);

% First pass
for i = 1 : K1
    id1 = i;
    v = parcels2(parcels1 == id1);
    [uniqs, counts] = countUniqueElements(v);
    [overlap, ind1] = max(counts);  
        
    id2 = uniqs(ind1);
    v = parcels1(parcels2 == id2);
    [uniqs, counts] = countUniqueElements(v);
    [~, ind2] = max(counts);
    cmp = uniqs(ind2);
    
    if cmp == id1    
        check1(id1) = 1;
        check2(id2) = 1;
        pairs(id1,1) = id1;
        pairs(id1,2) = id2;
        overlaps(id1) = overlap; 
        dices(id1) = (2 * overlap) / ...
                     (sum(parcels1 == id1) + sum(parcels2 == id2));
    end   
end

% Second pass
unassigned_ids = find(check2 == 0);
C = length(unassigned_ids);
for i = 1 : C
    id2 = unassigned_ids(i);
    v = parcels1(parcels2 == id2);
    [uniqs, counts] = countUniqueElements(v);
    [overlap, ind2] = max(counts);
    id1 = uniqs(ind2);
    
    if check1(id1) == 0  
        check1(id1) = 1;
        check2(id2) = 1;
        pairs(id1,1) = id1;
        pairs(id1,2) = id2;
        overlaps(id1) = overlap; 
        dices(id1) = (2 * overlap) / ...
                     (sum(parcels1 == id1) + sum(parcels2 == id2));
    end   
end

% Third pass
unassigned_ids = find(check1 == 0);
C = length(unassigned_ids);
for i = 1 : C
    id1 = unassigned_ids(i);
    v = parcels2(parcels1 == id1);
    [uniqs, counts] = countUniqueElements(v);
    [overlap, ind1] = max(counts);
    id2 = uniqs(ind1);
    
    if check2(id2) == 0  
        check1(id1) = 1;
        check2(id2) = 1;
        pairs(id1,1) = id1;
        pairs(id1,2) = id2;
        overlaps(id1) = overlap; 
        dices(id1) = (2 * overlap) / ...
                       (sum(parcels1 == id1) + sum(parcels2 == id2));
    end     
end


% Final pass
unassigned_ids = find(check2 == 0);
C = length(unassigned_ids);
for i = 1 : C
    id2 = unassigned_ids(i);
    v = parcels1(parcels2 == id2);
    [uniqs, counts] = countUniqueElements(v);
    [overlap, ind2] = max(counts);
    id1 = uniqs(ind2);
    
    if check1(id1) == 0  
        check1(id1) = 1;
        check2(id2) = 1;
        pairs(id1,1) = id1;
        pairs(id1,2) = id2;
        overlaps(id1) = overlap; 
        dices(id1) = (2 * overlap) / ...
                     (sum(parcels1 == id1) + sum(parcels2 == id2));
    end      
end

[sorted, idx] = sort(dices, 'Descend');
dices = sorted;
temp = pairs(:,1);
pairs(:,1) = temp(idx);
temp = pairs(:,2);
pairs(:,2) = temp(idx);

% Align matched clusters
id = 1;
for i = 1 : length(find(pairs(:,2) > 0)) 
    a = pairs(i,1);
    b = pairs(i,2);
    if b == 0
        continue;
    else
        Umatched(parcels1 == a) = id;
        Vmatched(parcels2 == b) = id;     
    end    
    id = id + 1;    
end

% Assign new ids to unmatched clusters
unassigned1 = setdiff(unique(parcels1), find(check1 == 1));
unassigned2 = setdiff(unique(parcels2), find(check2 == 1));
K1 = length(unassigned1);
K2 = length(unassigned2);

for i = 1 : max(K1,K2)
    if i <= K1
        Umatched(parcels1 == unassigned1(i)) = id;
    end
    if i <= K2
        Vmatched(parcels2 == unassigned2(i)) = id; 
    end
    id = id + 1;
end

% A final iteration to check if there is any remaining pair (very rare)
rems = find(dices == 0);
for i = 1 : length(rems)
    aa1 = Umatched == rems(i);
    aa2 = Vmatched == rems(i);
    dice = (2 * sum(aa1 == 1 & aa2 == 1)) / ...
           (sum(aa1 == 1) + sum(aa2 == 1));    
    if dice > 0
        dices(rems(i)) = dice;
    else
        break;
    end
end

dice = mean(dices);



