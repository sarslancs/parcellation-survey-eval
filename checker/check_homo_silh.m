root = '/vol/vipdata/data/HCP100/';
parrent = ['/vol/medic02/users/sa1013/Parcellation_Codes/' ... 
          '1_Two_Level_Clustering_Single_Subject/JOURNAL_Codes'];
      
addpath(genpath('/vol/medic02/users/sa1013/Parcellation_Codes/tools'));
addpath(genpath('/vol/medic02/users/sparisot/MATLAB/Parcellation'));
addpath(parrent);

dset = 2;
load(['subjectIDs100_set' num2str(dset) '.mat'], 'subjectIDs');
upto = 100;
fromSession = '1';

hems = ['L','R'];
lims = [29696 29716];
methods = cellstr([ 'GORDON';'ICA   ';'POWER ';'YEO   ';'BALD  ';'BROD  ';
                    'GLASS ';'COMPOS';'DEST  ';'DKT   ';'SHEN  ';'AAL   ';
                    'SLIC  ';'BLUM  ';'BELL  ';'WARD  ';'KMEANS';'NCUTS '; 
                    'GEO   ';'JOINT ';'GRASP ';'WARD  ';'KMEANS';'NCUTS ';
                    'GRAMPA';'BRAIN ']);
                
techs = cellstr([   'PRO   ';'PRO   ';'PRO   ';'PRO   ';'PRO   ';'PRO   ';
                    'PRO   ';'PRO   ';'MAJOR ';'MAJOR ';'PRO   ';'PRO   '; 
                    '2LEVEL';'2LEVEL';'2LEVEL';'2LEVEL';'2LEVEL';'2LEVEL';
                    'AVR   ';'ALL   ';'PCA   ';'PCA   ';'PCA   ';'PCA   '; 
                    'PCA   ';'PRO   ']);               
whichCorr = 'rS';
sim = 'dist';
maskThemOut = 1;
clossness   = 2;

HOMOs_ = zeros(length(Ks),length(methods));
HOMOs_New = zeros(length(Ks),length(methods));
SILHs = zeros(length(Ks),length(methods));
SILHsNew = zeros(length(Ks),length(methods));

for m = 1 : length(methods)
    if m == 6 || m == 8
        continue;
    end
    method = methods{m};
    tech =  techs{m};
    disp([method ' ' tech]);   
    loadName = fetch_group_name(method, tech, 1, hem);
    load(loadName, 'PARCELs'); 
    Ks = PARCELs.Ks;
    tic
    for k = 1 : length(Ks)  
        parcellation = PARCELs.parcels{k};
        neigh = PARCELs.neighSet{k};
        [ ~, homos_, ~] = silhouette_analysis_craddock( parcellation, Z, neigh, vertices, 2);
        HOMOs_(k,m) = mean(homos_);
        [ silhs, ~, ~ ] = silhouette_analysis_yeo( parcellation, corrs, neigh, vertices, sim, clossness, maskThemOut);         
        SILHs(k,m) = mean(silhs);
        
        [ homogeneities ] = homogeneity( parcellation, Z );
        HOMOs_New(k,m) = mean(homogeneities);
        [ silhouettes ] = silhouette_coef( parcellation, 1-corrs, neigh );
        SILHsNew(k,m) = mean(silhouettes);
    end
    toc 
    clear PARCELs 
end
  
