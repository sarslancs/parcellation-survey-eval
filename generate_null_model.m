function [ nullpar, moved ]= generate_null_model( parcels, sphere, mask )
%GENERATE_NULL_MODEL Generate a null model (parcellations).
%   For each hemisphere, a parcellation is projected onto the standard
%   spherical surface (provided by the HCP) and each point in this sphere  
%   is randomly rotated around the x, y, and z axes. This process moves  
%   each parcel to a new location on the cortical surface without altering 
%   their relative positions. Such null models can then be used for 
%   homogeneity analysis. 
%
%   INPUT
%   =====
%   parcels: A parcellation.
%   sphere: A spherical surface model (can be obtained from a 
%           sphere.32k_fs_LR.surf.gii file)
%   mask: A binary cortical mask, in which a vertex v = 0 if v is in 
%         medial wall, othwerwise v = 1 (can be acquired from an 
%         atlasroi.32k_fs_LR.shape.gii file)
%
%   OUTPUT
%   ======
%   nullpar: A randomly rotated parcellation.
%   moved: Labels of parcels that were moved to the medial wall, and thus,
%            will be discarded from analysis.
%
%   USAGE
%   =====
%   [ NULLPAR, MOVED ] = GENERATE_NULL_MODEL( PARCELS, SPHERE, MASK ) 
%   returns a N-by-1 null parcellation (NULLPAR), and a list of parcel 
%   indices (MOVED) that are within the medial wall after rotation.
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


thetas = randi([30, 150],1,3); % The range of rotations

rotated = ((sphere.vertices * rotx(thetas(1))) * ...
                              roty(thetas(2))) * ...
                              rotz(thetas(3));
IDX = knnsearch(sphere.vertices,rotated);
mm = zeros(size(mask));
mm(logical(mask)) = parcels;
nullpar = mm(IDX);
temp = nullpar;
nullpar(~mask) = [];
temp(mask == 1) = 0;
moved = unique(nonzeros(temp)); % Do not compute homogeneity for the 
% parcels that end up in the medial wall. Instead assign these parcels 
% the average homogeneity of all other parcels that were rotated into 
% valid cortical areas, as suggested in Gordon et al. (2016), Cereb Cortex


