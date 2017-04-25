% Continue later!

% Subject ID (first subject in Dataset 1)
subjectID = '100307' ;

% Parameters 

hem = 'L'; % Hemisphere
session = '1'; % Session
method = 'Arslan'; % Parcellation method

loadFrom = ['Subject/' method '/']; 
loadName = [loadFrom subjectID '_' method '_' session '_' hem];     
load(loadName, 'parcels');

% Get a parcellation with the desired resolution
parcellation = parcels(:,19); % i.e. K = 100
           
% Load subject's myelin map from the HCP complemntary files in 
% MNINonLinear/fsaverage_LR32k/
thresh = 1;
session = '1';

        disp(subjectID); 
        load([root subjectID '/processed/' subjectID '_atlasroi_cdata_' hem], 'cdata');
        myelFrom = [root subjectID '/structural/MNINonLinear/fsaverage_LR32k/'];
        g = gifti([myelFrom subjectID '.' hem '.MyelinMap.32k_fs_LR.func.gii']);
        myel = g.cdata;
        myelin = myel(cdata>0);
        
%         file:///vol/medic02/users/sparisot/ClionProjects/openGm_maxflow/build/100307/L/myelin_100307_25_initR10.txt

        % Compute the average map
        myelFrom = ['/vol/medic02/users/sparisot/ClionProjects/openGm_maxflow/build/' subjectID  ...
                     '/' hem '/myelin_' subjectID '_25_initR10.txt'];
        myel = dlmread(myelFrom);
        avM = myel';
        avM = avM(cdata>0); 
        
        relabeled = relabelParcels(avM,max(avM));
        if thresh
            for j = 1 : max(relabeled)
                means = mean(myelin(relabeled == j));
                disp(means)
                if means < median(myelin)
                    relabeled(relabeled == j) = 0;
                end 
            end
        end
            
        for m = 1 : length(methods)
            met = methods{m};
            disp(met);        
            loadName =  fetch_single_name(met, subjectID, params, session, hem, root);
            try
                load(loadName, 'PARCELs'); 
            catch
                error(['Something is wrong! Better check subject ' ...
                num2str(i), '-' met ', for session ' num2str(ii)', ...
                hem ' hem '.']);
            end

            tic
            for k = 1 : length(Ks)  
                if ~isempty(PARCELs.parcelSet{k})
                    parcellation = PARCELs.parcelSet{k};
                    if size(parcellation,1) == 1
                        parcellation = parcellation';
                    end
                    [ dices, couples, lab1, lab2 ] = diceOverlapJoined( parcellation, relabeled);
                    MYELsDICE(i,k,m) = mean(dices);
                else
                    MYELsDICE(i,k,m) = NaN; 
                end
            end
            toc 
            clear PARCELs 
        end
        save([parrent '/results/SINGLE_MYELIN_DICE_set_' num2str(dset) '_' whichMyel '_thresh_' num2str(thresh) '_ALL_' hem ], 'MYELsDICE');  
    end
end