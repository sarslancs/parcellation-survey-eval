function [ uniqs, counts ] = count_unique_elements( V )
%COUNT_UNIQUE_ELEMENTS Counts the number of unique elements in the array
%   Given an array of v counts the number of unique elements. 
%   
%   INPUT
%   =====
%   V: A Nx1 column vector, where N is the number of elements. 
%
%   OUTPUT
%   ======
%   uniqs: M x 1 vector, where M is the number of unique elements in 
%          the input vector. 
%   counts: M x 1 vector, in which counts_i corresponds to the number of 
%   appearances (i.e. counts) of the value indexed at i in uniqs. 
%   
%   CAUTION
%   =======
%   V must be a column vector, otherwise the following error may be 
%   encountered. counts is descendingly sorted.
%   'Error using vertcat
%   Dimensions of matrices being concatenated are not consistent.'
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

    y = sort(V);
    p = find([true;diff(y)~=0;true]);
    uniqs = y(p(1:end-1));
    counts = diff(p);
    [counts, idx] = sort(counts, 'descend');
    uniqs = uniqs(idx);   
end

