# Brain Parcellation Survey: Evaluation Code
This repo holds the evaluation code (mainly written in Matlab) used in the brain parcellation survey "Human brain mapping: 
A systematic comparison of parcellation methods for the human cerebral cortex", in which a large-scale and systematic 
comparison is carried out between the state-of-the-art connectivity-driven, anatomical, and random parcellation methods. 

Please cite the corresponding paper (if using the code) as: 

Arslan, S., Ktena, S.I., Makropoulos, A., Robinson, E.C., Rueckert, D., Parisot, S., "Human brain mapping: A systematic 
comparison of parcellation methods for the human cerebral cortex", NeuroImage (2017), http://dx.doi.org/10.1016/j.neuroimage.2017.04.014

## Some useful links you may want to see _in prior_ using the codes

Survey page: https://biomedia.doc.ic.ac.uk/brain-parcellation-survey/

Download parcellations from: https://imperialcollegelondon.box.com/s/g5q0kyvpqdha5jgofhmiov9ws1ao0hi0

Preprint of manuscript: https://www.doc.ic.ac.uk/~sa1013/pub/2017-neuroimage-brain-parcellation-survey.pdf

Data reference manual V2: https://biomedia.doc.ic.ac.uk/wp-content/uploads/sites/95/2017/03/Reference-Manual.pdf 

## Usage

Simply follow in-file instructions. Reading the Data and Experimental Setup sections in the paper may also 
help understand the code.

## Data
Generation and evaluation of parcellations are carried out using the publicly available dataset from the Human Connectome 
Project (HCP). Some functions may be dependednt on the HCP data structures and require you to download certain files. 
However, each script can also be adapted to any other dataset, by following the steps throughout the codes. If you would like 
to use the HCP data, it can be downloaded from https://db.humanconnectome.org 

## Parcellations
To follow the experimental pipeline described in the paper, all parcellations are made publicly available at
https://biomedia.doc.ic.ac.uk/brain-parcellation-survey/ Please see the data reference manual for a quick walkthrough.


