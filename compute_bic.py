"""
Created Apr 2017

@author: Salim Arslan (name.surname@imperial.ac.uk)
"""

"""
Given task fMRI contrast maps and parcellations, computes the BIC scores for 
the evaluation of parcellation accuracy with respect to task activation. An 
example processing pipeline is provided at the end.

The analysis code used in this script is described in "Which fMRI clustering 
gives good brain parcellations?" (doi.org/10.3389/fnins.2014.00167) and 
available at https://github.com/bthirion/frontiers_2014. We have included it 
to our evaluation pipelines in the brain parcellation survey, "Human brain 
mapping: A systematic comparison of parcellation methods for the human cerebral 
cortex" (doi.org/10.1016/j.neuroimage.2017.04.014) 

For the parcellation data and reference manual visit the survey page: 
https://biomedia.doc.ic.ac.uk/brain-parcellation-survey/ 

 
"""

# Import necessary modules
from mixed_effects_parcel import parameter_map 

# This is just a simple function that calls the original BIC code     
def compute_bic(X, parcels):
    """Return the BIC of a parcellation w.r.t. a contrast map

    Parameters
    ==========
    X: numpy array of shape(n_vertices, n_contrasts, n_subjects) contrast map 
       matrix for a set of subjects (e.g. all subjects in Dataset 2) and 
       contrasts, which can be activation maps (e.g. z scores) from a subject-
       level task fMRI analysis.  
    parcels: numpy array of shape (n_vertices) an index array describing the 
             parcellation, where label index starts from 0 (if the provided
             survey parcellations are used, labels must be decremented by 1).

    Returns
    =======
    bic: float, BIC score
    """        
    bic = 0
    for i in range(X.shape[1]): #For each contrast, compute BIC and sum over
        ll_, mu_, sigma1_, sigma2_, bic_ = parameter_map(X[:,i], 
                                                         parcels, 
                                                         null=False)
        bic += bic_.sum()
        
    return bic
    
    
""" An example analysis pipeline

import numpy as np

# Parameters
hem = 'L'               # Left hemisphere
n_subjects = 100        # All subjects in Dataset 2
n_voxels = 29696        # Number of vertices in L hemisphere
n_contrasts = 86        # Number of contrasts
                                 
# Initialise contrast map matrix
X = np.zeros((n_voxels, n_contrasts, n_subjects)) # Contrasts of all subjects 
are stored in X, a numpy array of shape(n_vertices, n_contrasts, n_subjects)    

# Populate contrast matrix X
for contrast in contrasts: # For each contrast (e.g. MOTOR-LF-AVG)    
    # Populate contrast matrix X

# Run BIC analysis for each method/resolution        
for method in methods: # For each method
    for K in enumerate(Ks):  # For each resolution 
        # Get parcels and run compute_bic 
        bic = compute_bic(X, parcels)      

# Tip: In order to use the provided parcellations downloaded from 
# https://biomedia.doc.ic.ac.uk/brain-parcellation-survey/, see 
# load_all_group_parcellations.py or load_all_single_subject_parcellations.py,                   
# located in the "Scripts" folder

"""