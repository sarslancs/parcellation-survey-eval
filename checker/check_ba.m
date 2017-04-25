parrent = ['/vol/medic02/users/sa1013/Parcellation_Codes/' ... 
          '1_Two_Level_Clustering_Single_Subject/JOURNAL_Codes'];
      
addpath(genpath('/vol/medic02/users/sa1013/Parcellation_Codes/tools'));

root = '/vol/vipdata/data/HCP100/';

dset = 2;
load(['subjectIDs100_set' num2str(dset) '.mat'], 'subjectIDs');

hems = ['L','R'];
nMesh = 32492;

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
               
Ks = [10:5:100 125:25:250];

for p = 1 : 2
    hem = hems(p);

    OVERLAPs = zeros(length(methods), length(Ks));
    OVERLAPsNew = zeros(length(methods), length(Ks));

    root = '/vol/vipdata/data/HCP100/';
    subjectID = '100307';
    getFrom = [root subjectID '/processed/'];        
    load([getFrom subjectID '_atlasroi_cdata_' hem], 'cdata');
    
    singles = zeros(32492,length(subjectIDs));
    for i = 1 : length(subjectIDs)
        subjectID = num2str(subjectIDs(i));
        if dset == 2
            if i <= 50
                root = '/vol/dhcp-hcp-data/twins_data/';
            else
                root = '/vol/dhcp-hcp-data/HCP50M/';
            end
        end

        disp(subjectID); 
        
        getFrom = [root subjectID '/processed/'];        
        load([getFrom subjectID '_atlasroi_cdata_' hem], 'cdata');

        %%% Brodmann individual subject
        from = ([root subjectID '/structural/MNINonLinear/fsaverage_LR32k/']);
        bas = gifti([from subjectID '.' hem '.BA.32k_fs_LR.label.gii' ]);
        tmp=bas.cdata(:,1);
        
        % 1 -> BA1
        % 2 -> BA2
        % 3 -> BA3a
        % 4 -> BA3b
        % 5 -> BA4a
        % 6 -> BA4p
        % 7 -> BA6
        % 8 -> BA44
        % 9 -> BA45
        % 10 -> V1
        % 11 -> V2
        % 12 -> MT
        % 13 -> BA35,36
        
        tmp(tmp==2 | tmp==3 | tmp==4) = 1; % -> primary somatosensory cortex 
        tmp(tmp==6) = 5; %->  motor cortex (saved as BA7)
        tmp(tmp==11) = 10; %-> visual cortex
        
%         tmp(tmp==6) = 5; %-> primary motor cortex
%         tmp(tmp==4 | tmp==3 | tmp==5 | tmp==6)=1;
%         tmp(tmp==11)=10;
%         tmp=tmp+1;
        singles(:,i) = tmp;
    end
    atlas = majority_voting_on_sets( singles, nMesh);
    atlas = atlas(cdata > 0);
    ids = unique(atlas);
       
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
            if ~isempty(parcellation)

                [ dices, ids, lab1, lab2 ] = diceOverlapJoined( parcellation, atlas);
                OVERLAPs(m,k) = mean(dices);

                load(['atlas/atlas_group_avr_' hem '.mat'], 'ba')
                [ dice, ~ ] = dice_atlas_align( parcellation, ba);
                OVERLAPsNew(m,k) = dice;

            else
                OVERLAPs(m,k) = NaN;
            end
        end
    end
%     save([parrent '/results/GROUP_OVERLAPs_BA8_ALL_' hem ], 'OVERLAPs')
end