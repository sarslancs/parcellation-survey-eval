function [ dice, joined ] = dice_atlas_align( parcels, atlas )
%DICE_ATLAS_ALIGN Alignment (similarity) of a parcellation with an atlas.
%   Computes a Dice score to assess the similarity of a parcellation 
%   to an atlas (e.g. Brodmann areas or any other reference parcellation).
%   Similar to dice_coef_joined, an iterative procedure is performed to 
%   indentify matching parcels before computing pairwise Dice overlaps. 
%   Parcels matching with the same atlas region are joined (merged), before 
%   computing pairwise Dice overlaps. A global Dice score is obtained via 
%   averaging.
%
%   INPUT
%   =====
%   parcels: Input parcellation
%   atlas: Reference parcellation (e.g. Brodmann atlas)
%
%   OUTPUT
%   ======
%   dice: Average Dice score computed between parcels and atlas
%   joined: Input parcels with new labels after matching and joining 
%   
%   USAGE
%   =====
%   [ DICE, JOINED] = dice_atlas_align( PARCELS, ATLAS) returns the DICE
%   score that measures the similarity between the target and reference
%   parcellations (PARCELS and ATLAS, respectively). PARCELS can be a 
%   subject-level or groupwise parcellation. ATLAS can be the Brodmann 
%   areas or any other reference parcellation.
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
%
%   
%   TUTORIAL
%   ========
%   For example, to load the group avereage reference atlases for the left 
%   hemisphere, run: 
%   
%   load('atlas/atlas_group_avr_L')  
%   
%   This will return three parcellations, namely MYELIN, MYELIN_THR and BA. 
%   MYELIN is the original coarse-resolution myelin parcellation described  
%   in the survey paper (see Fig. 3b). MYELIN_THR is the thresholded myelin
%   parcellation, in which only areas with high myelination are retained
%   (this is the reference parcellation used for measuring alignment to 
%   structured patterns of myelination). BA is a parcellation with labels
%   of several Brodmann areas (as shown in Fig. 3a). The corresponding
%   labels in Fig. 3a are mapped to BA as follows:
%
%   BA[3,1,2]   -> 1
%   BA4         -> 2
%   BA6         -> 3
%   BA44        -> 4
%   BA45        -> 5
%   BA17        -> 6
%   MT          -> 7
%   BA[35,36]   -> 8
%

parcelsIN = parcels;
parcels(atlas == 0) = 0;

newparcels = zeros(size(parcels));
joined = zeros(size(atlas));

clusters1 = unique(parcels);
clusters1(clusters1 == 0) = [];
K1 = length(clusters1);

clusters2 = unique(atlas);
clusters2(clusters2 == 0) = [];
K2 = length(clusters2);

ind = 1;
for i = 1 : length(clusters1);
   id1 = clusters1(i);
   newparcels(parcelsIN == id1) = ind;
   ind = ind + 1;
end

ind = 1;
for i = 1 : length(clusters2);
   id1 = clusters2(i);
   joined(atlas == id1) = ind;
   ind = ind + 1;
end

parcels = newparcels;
atlas = joined;

joined = zeros(size(atlas));

check1 = zeros(K1, 1);
check2 = zeros(K2, 1);

pairs = zeros(K2,2);
dices = zeros(K2,1);

% Join them all
for i = 1 : length(clusters2);
    v = parcels(atlas == i);
    [ instances, values ] = countUniqueElements(v);
    sums = arrayfun(@(x) sum(parcels == x), instances);
    res = (values./sums) >= .5;    
    if sum(res) > 1
        merged = instances(res);
        newid = merged(1);        
        for m = 1 : length(merged)
            parcels(parcels == merged(m)) = newid;
        end
    end
end


clusters2 = nonzeros(unique(atlas));
K2 = length(clusters2);

% Match them all
for i = 1 : K2
    id1 = clusters2(i);
    v = parcels(atlas == id1);
    [uniqs, counts] = countUniqueElements(v);
    sums = arrayfun(@(x) sum(parcels == x), uniqs);
    [~, ind1] = max(counts./sums); 
        
    id2 = uniqs(ind1);
    if id2 > 0
        v = atlas(parcels == id2);
        [uniqs, counts] = countUniqueElements(v);
        sums = arrayfun(@(x) sum(atlas == x), uniqs);
        [~, ind2] = max(counts./sums);
        cmp = uniqs(ind2);
        if cmp == id1    
            check1(id2) = 1;
            check2(id1) = 1;
            pairs(id1,1) = id1;
            pairs(id1,2) = id2;
            dices(id1) = (2 * sum(atlas == id1 & parcels == id2)) / ...
                         (sum(atlas == id1) + sum(parcels == id2));
        end
    end      
end

% Check if there are still parcels to match (very rare)
unassigned_ids = nonzeros(setdiff(unique(atlas), find(check2 == 1)));
K1 = length(unassigned_ids);
for i = 1 : K1
    id2 = unassigned_ids(i);
    v = parcels(atlas == id2);
    [uniqs, counts] = countUniqueElements(v);
    sums = arrayfun(@(x) sum(parcels == x), uniqs);
    [~, ind2] = max(counts./sums); 
    id1 = uniqs(ind2);
    if check1(id1) == 0  
        check1(id1) = 1;
        check2(id2) = 1;
        pairs(id2,1) = id2;
        pairs(id2,2) = id1;
        dices(id2) = (2 * sum(atlas == id2 & parcels == id1)) / ...
                     (sum(parcels == id1) + sum(atlas == id2));
    end      
end

for i = 1 : length(pairs)
    if pairs(i,1) ~= 0
        joined(parcels == pairs(i,2)) = pairs(i,1);
    end
end

dice = mean(dices);




