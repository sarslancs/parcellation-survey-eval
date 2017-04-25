function [ dice, Ujoined, Vjoined ] = dice_coef_joined( U, V )
%DICE_COEF_JOINED Joined Dice coefficient.
%   Computes the Dice score to assess the similarity of parcellations in a
%   scan-to-scan or group-to-group setting. An iterative procedure is 
%   performed to indentify matching parcels before computing pairwise 
%   Dice overlaps. Overlapping parcels are joined to reduce impact of 
%   over-parcellation. A global Dice score is obtained via averaging.
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
%   Ujoined: Parcellation U, with new labels after matching and joining 
%   Vjoined: Parcellation V, with new labels after matching and joining 
%
%   USAGE
%   =====
%   DICE = DICE_COEF_JOINED( U, V ) returns a DICE value, which indicates  
%   the similarity between parcellations U and V. U and V can be (1) 
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

Ujoined = zeros(size(parcels1));
Vjoined = zeros(size(parcels2));

check1 = zeros(K1, 1);
check2 = zeros(K2, 1);

overlaps = zeros(max(K1,K2),1);
pairs = zeros(max(K1,K2),2);
dices = zeros(max(K1,K2),1);

% JOIN: First pass
for i = 1 : K1;
    v = parcels2(parcels1 == i);
    [ instances, values ] = countUniqueElements(v);
    sums = arrayfun(@(x) sum(parcels2 == x), instances);
    res = (values./sums) >= .5;
    
    if sum(res) > 1
        merged = instances(res);
        newid = merged(1);           
        for m = 1 : length(merged)
            parcels2(parcels2 == merged(m)) = newid;
        end
    end
end

% JOIN: Second pass
clusters2 = unique(parcels2);
for i = 1 : length(clusters2);
    v = parcels1(parcels2 == clusters2(i));
    [ instances, values ] = countUniqueElements(v);
    sums = arrayfun(@(x) sum(parcels1 == x), instances);
    res = (values./sums) >= .5;
    
    if sum(res) > 1
        merged = instances(res);
        newid = merged(1);           
        for m = 1 : length(merged)
            parcels1(parcels1 == merged(m)) = newid;
        end
    end
end


% MATCH: First pass
clusters1 = unique(parcels1);
for i = 1 : length(clusters1)
    id1 = clusters1(i);
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

% MATCH: Second pass
unassigned_ids = setdiff(unique(parcels2), find(check2 == 1));
for i = 1 : length(unassigned_ids)
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

% MATCH: Third pass
unassigned_ids = setdiff(unique(parcels1), find(check1 == 1));
for i = 1 : length(unassigned_ids)
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


ids = find(pairs(:,1) == 0 & pairs(:,2) == 0);
pairs(ids,:) = [];
dices(ids, :) = [];

% Align matched clusters
[sorted, idx] = sort(dices, 'Descend');
dices = sorted;
temp = pairs(:,1);
pairs(:,1) = temp(idx);
temp = pairs(:,2);
pairs(:,2) = temp(idx);

id = 1;
for i = 1 : length(pairs(:,2)) 
    a = pairs(i,1);
    b = pairs(i,2);
    if b == 0
        continue;
    else
        Ujoined(parcels1 == a) = id;
        Vjoined(parcels2 == b) = id;     
    end    
    id = id + 1;    
end

% Assign labels to unmatched parcels
unassigned1 = setdiff(unique(parcels1), find(check1 == 1));
unassigned2 = setdiff(unique(parcels2), find(check2 == 1));
C1 = length(unassigned1);
C2 = length(unassigned2);

for i = 1 : max(C1,C2)
    if i <= C1
        Ujoined(parcels1 == unassigned1(i)) = id;
    end
    if i <= C2
        Vjoined(parcels2 == unassigned2(i)) = id; 
    end
    id = id + 1;
end

% Account for the unmatched parcels in the global Dice score
dices(end + 1:max(max(Vjoined), max(Ujoined))) = 0;
dice = mean(dices);



