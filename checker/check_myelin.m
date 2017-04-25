parrent = ['/vol/medic02/users/sa1013/Parcellation_Codes/' ... 
          '1_Two_Level_Clustering_Single_Subject/JOURNAL_Codes'];
      
addpath(genpath('/vol/medic02/users/sa1013/Parcellation_Codes/tools'));
addpath(genpath('/vol/medic02/users/sparisot/MATLAB/Parcellation'));
addpath(parrent);

dset            = 2;
nSub            = 100;
hems            = ['L','R'];
lims            = [29696,29716];
Ks              = [10:5:100 125:25:250];
whichMyel       = 'MyelinMap_BC'; %-> Confirm with Sarah
thresh          = 1;

root = '/vol/vipdata/data/HCP100/';
subjectID = '100307';
getFrom = [root subjectID '/processed/'];

load(['subjectIDs100_set' num2str(dset) '.mat'], 'subjectIDs');
            
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

                
for p = 1 : 1  
    MYELs = zeros(length(methods), length(Ks) ); 
    MYELsNew = zeros(length(methods), length(Ks) ); 
    hem = hems(p);   
    disp([hem ' - sim: ' whichMyel])
    
    % Compute the average map
    avM = zeros(lims(p),1);
    for i = 1 : nSub
        subjectID = num2str(subjectIDs(i));
        if dset == 2
            if i <= 50
                root = '/vol/dhcp-hcp-data/twins_data/';
            else
                root = '/vol/dhcp-hcp-data/HCP50M/';
            end
        end

        disp(subjectID); 
        
        load([root subjectID '/processed/' subjectID '_atlasroi_cdata_' hem], 'cdata');
        myelFrom = [root subjectID '/structural/MNINonLinear/fsaverage_LR32k/'];
        g = gifti([myelFrom subjectID '.' hem '.' whichMyel '.32k_fs_LR.func.gii']);
        myel = g.cdata;
        avM = avM + myel(cdata>0); 
    end
    myelin = avM / nSub;
    
    
    load([root subjectID '/processed/' subjectID '_atlasroi_cdata_' hem], 'cdata');
    myelFrom = ['/vol/medic02/users/sparisot/ClionProjects/openGm_maxflow/build/GrAverage100PCA' ...
                 num2str(dset) '/' hem '/myelin_GrAverage100PCA' num2str(dset) '__'];
    myel = dlmread(myelFrom);
    avM = myel';
    avM(cdata == 0) = [];
    
    if thresh
        for i = 1 : max(avM)
            means = mean(myelin(avM == i));
            if means < median(myelin)
                avM(avM == i) = 0;
            end
        end
    end
            
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
        for k = 1 : length(Ks)
            parcellation = PARCELs.parcels{k};
            [ dices, ids, lab1, ~ ] = diceOverlapJoined( parcellation, avM);
            
            MYELs(m,k) = mean(dices);
            
            load(['atlas/atlas_group_avr_' hem '.mat'], 'myelin_thr')
            [ dice, ~ ] = dice_atlas_align( parcellation, myelin_thr);
            MYELsNew(m,k) = dice;
        end
        clear PARCELs        
    end  
%     save([parrent '/results/GROUP_MYELIN_DICE_set_' num2str(dset) '_' whichMyel '_thresh_' num2str(thresh) '_ALL_' hem ], 'MYELs');  
end