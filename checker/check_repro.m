parrent = ['/vol/medic02/users/sa1013/Parcellation_Codes/' ... 
          '1_Two_Level_Clustering_Single_Subject/JOURNAL_Codes'];
addpath(genpath('/vol/medic02/users/sa1013/Parcellation_Codes/tools'));
addpath(parrent);

hems = ['L','R'];

methods = cellstr([ 'SLIC  ';'BLUM  ';'BELL  ';'WARD  ';'KMEANS';'NCUTS '; 
                    'GEO   ';'JOINT ';'GRASP ';'WARD  ';'KMEANS';'NCUTS ';
                    'GRAMPA']);
                
techs = cellstr([   '2LEVEL';'2LEVEL';'2LEVEL';'2LEVEL';'2LEVEL';'2LEVEL';
                    'AVR   ';'ALL   ';'PCA   ';'PCA   ';'PCA   ';'PCA   '; 
                    'PCA   ']);
                
Ks = [10:5:100 125:25:250];

for p = 1 : 1
    DICEs = zeros(length(methods), length(Ks));
    DICEJs = zeros(length(methods), length(Ks));
    DICEsNew = zeros(length(methods), length(Ks));
    DICEJsNew = zeros(length(methods), length(Ks));
    hem = hems(p);   
    disp(hem);

    for m = 1 : length(methods)
        method = methods{m};
        tech = techs{m};
        disp([method ' ' tech]);   

        loadName = fetch_group_name(method, tech, 1, hem);

        load(loadName, 'PARCELs'); 
        PARCELs1 = PARCELs; clear PARCELs
        Kin = PARCELs1.Ks;

        loadName =  fetch_group_name(method, tech, 2, hem);
        load(loadName, 'PARCELs'); 
        PARCELs2 = PARCELs; clear PARCELs 

        tic
        for k = 1 : length(Kin)  
            [ dices, ~, ~ ] = alignLabels3( PARCELs1.parcels{k}, PARCELs2.parcels{k}, 1 ); 
            DICEs(m,k) = mean(dices); 
            [ dices, ~, ~ ] = alignLabelsJoined3( PARCELs1.parcels{k}, PARCELs2.parcels{k}, 1 ); 
            DICEJs(m,k) = mean(dices); 
            
            DICEsNew(m,k) = dice_coef( PARCELs1.parcels{k}, PARCELs2.parcels{k}); 
            DICEJsNew(m,k) = dice_coef_joined( PARCELs1.parcels{k}, PARCELs2.parcels{k}); 
            
        end
        toc 

        clear PARCELs1 PARCELs2 
    end
%     save([parrent '/results/GROUP_DICEs_ALL_' hem ], 'DICEs')
%     save([parrent '/results/GROUP_DICEJs_ALL_' hem ], 'DICEJs')
end
    
    
        