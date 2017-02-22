%ACE_Scanning_Analyzer_v1_0.m
%%  Main script to run an ACE Scanning microarray analysis from. 
%   See to the readme file for help and first steps.


%% Clear the matalb environment
clc;        % Clear the command window
clear all;  % Clear workspace
close all; % Close all figures
workspace;  % Make sure the workspace panel is showing




%% 1. SET THE DIRECTORY - update the following line to point to the location of this folder on your system
cd('/Users/jmunzar/GitHub/ACE-Scanning/');

% Add required code to Matlab search path
cd 'Code/'
addpath(genpath(cd)); 
cd ..





%% 2. SELECT THE ACE SCANNING DATASET TO ANALYZE
% Uncomment a single dataset and run the script to run the analysis


%%% ATP DNA 1 (ATP DNA - 10^12 density - 2015.08.04 - ATP)
cd ('RawDataFiles/ATP_DNA_1_ATP_12density'); load('ATP_DNA_1_ATP_12density.mat'); cd ../..;


%%% ATP DNA 2 (ATP DNA - Negative Control - 2016.07.13 - GTP)
% cd ('RawDataFiles/ATP_DNA_2_GTP'); load('ATP_DNA_2_GTP.mat'); cd ../..;


%%% ATP DNA 3 (ATP DNA - 10^11 density - 2016.02.25 - ATP)
% cd ('RawDataFiles/ATP_DNA_3_11density'); load('ATP_DNA_3_11density.mat'); cd ../..;


%%% ATP DNA 4 (ATP DNA - 10^10 density - 2016.02.25 - ATP)
% cd ('RawDataFiles/ATP_DNA_4_10density'); load('ATP_DNA_4_10density.mat'); cd ../..;


%%% ATP DNA 5 (ATP DNA - 10^9 density - 2016.02.25 - ATP)
% cd ('RawDataFiles/ATP_DNA_5_9density'); load('ATP_DNA_5_9density.mat'); cd ../..;


%%% ATP DNA 6 (ATP DNA - Quantitaitive Slide 1 - 2016.07.06 - ATP)
% cd ('RawDataFiles/ATP_DNA_6_Quant1'); load('ATP_DNA_6_Quant1.mat'); cd ../..;


%%% ATP DNA 7 (ATP DNA - Quantitaitive Slide 2 - 2016.07.06 - ATP)
% cd ('RawDataFiles/ATP_DNA_7_Quant2'); load('ATP_DNA_7_Quant2.mat'); cd ../..;


%%% ATP RNA 1 (ATP RNA - 10^12 density - 2016.02.28 - ATP)
% cd ('RawDataFiles/ATP_RNA_1'); load('ATP_RNA_1.mat'); cd ../..;


%%% Cocaine DNA 1 (Cocaine DNA - 10^12 density - 2016.02.06 - Cocaine)
% cd ('RawDataFiles/Cocaine_DNA_1'); load('Cocaine_DNA_1.mat'); cd ../..;


%%% Thrombin 1 Green (Thrombin DNA - 10^10 density - 2016.02.28 - Thrombin in Sodium Buffer)
% cd ('RawDataFiles/Thrombin_1_Sodium_Green'); load('Thrombin_1_Sodium_Green.mat'); cd ../..;


%%% Thrombin 1 Red (Thrombin DNA - 10^10 density - 2016.02.28 - Thrombin in Sodium Buffer)
% cd ('RawDataFiles/Thrombin_1_Sodium_Red'); load('Thrombin_1_Sodium_Red.mat'); cd ../..;


%%% Thrombin 2 Green (Thrombin DNA - 10^10 density - 2016.07.13 - Thrombin in Potassium Buffer)
% cd ('RawDataFiles/Thrombin_2_Potassium_Green'); load('Thrombin_2_Potassium_Green.mat'); cd ../..;


%%% Thrombin 2 Red (Thrombin DNA - 10^10 density - 2016.07.13 - Thrombin in Potassium Buffer)
% cd ('RawDataFiles/Thrombin_2_Potassium_Red'); load('Thrombin_2_Potassium_Red.mat'); cd ../..;


%%% add Riboswitch 1 (add RNA - 10^10 density - 2016.07.07 - Adenine and 4mM Mg at 23C)
% cd ('RawDataFiles/addRiboswitch_1_23C_4mM'); load('addRiboswitch_1_23C_4mM.mat'); cd ../..;


%%% add Riboswitch 2 (add RNA - 10^10 density - 2016.08.16 - Adenine and 12 mM Mg at 23C)
% cd ('RawDataFiles/addRiboswitch_2_23C_12mM'); load('addRiboswitch_2_23C_12mM.mat'); cd ../..;


%%% add Riboswitch 3 (add RNA - 10^10 density - 2016.08.11 - Adenine with 0 mM Mg at 10C)
% cd ('RawDataFiles/addRiboswitch_3_10C_0mM'); load('addRiboswitch_3_10C_0mM.mat'); cd ../..;


%%% add Riboswitch 4 (add RNA - 10^10 density - 2016.08.11 - Adenine with 12 mM Mg at 10C)
% cd ('RawDataFiles/addRiboswitch_4_10C_12mM'); load('addRiboswitch_4_10C_12mM.mat'); cd ../..;










%% Create Results Folder if not already made

cd('CompiledResults/')
mkdir(workingfolder); % Make a folder to save the results of the analysis
cd ..



%% Change FLAGS for analysis (Defaults were used for publication)

%%% rmOutliersBKGND FLAG:
%   1:  Remove outliers based on background data 
%       (stdev of mean background signal)
%   0:  Do not remove any outliers based on background
FLAG_rmOutliersBKGND = 1 %Default = 1
FLAG_rmOutliersBKGND_sigma = 3 %Default = 3


%%% rmOutliersMorphology FLAG:
%   1:  Remove outliers based on spot morphology
%       (stdev of spot signal)
%   0:  Do not remove any outliers based on morphology
FLAG_rmOutliersMorphology = 1 %Default = 1
FLAG_rmOutliersMorphology_sigma = 3 %Default = 3


%%% rmLowHybdata FLAG:
%%% 2:  Remove poorly hybridized ACEs ( Hyb signal < hybSignalFloor 
%%%     and/or Assay signal < assaySignalFloor )
%%% 1:  Remove poorly hybridized Hyb ACEs only (Hyb signal < hybSignalFloor)
%%% 0:  Use all ACEs
FLAG_rmLowHybdata = 2 %Default = 2
FLAG_rmLowHybdata_HybFloor = 200 %Default = 200
FLAG_rmLowHybdata_AssayFloor = 100 %Default = 100






%%  IMPORT MICROARRAY DATA
%%%
%%%

cd('RawDataFiles/');
cd(workingfolder);

% Get spot location, sorting information from the gal file:
[N_Rows, N_Columns, N_SubArrays, N_ACEs, N_ReplicatesPerACE, ...
    SpotIDnumber_Left, SpotIDnumber_Right, ...
    sort_indices_Left, sort_indices_Right] ...   
    = extractGalFile(filename_Gal, cd);

% Import raw data into structures:
[Hyb_rawdata, Assay_rawdata] = extractArrayData(filename_Hyb, ...
    filename_Assay, cd, N_SubArrays);

cd ../..




%% ORIENT DATA IN MICROARRAY LAYOUT FOR DATA VIEWING
% Use OrientedData to look at raw microarray data

OrientedData(1:N_Rows,1:N_Columns,1:N_SubArrays) = NaN;
counter = 0;
for i = 1:N_Rows
    for j = 1:N_Columns
        counter = counter+1;
        for k = 1:N_SubArrays
            OrientedData(i,j,k) = ...
                    Hyb_rawdata.stdevSignal(counter,k) / ...
                    Hyb_rawdata.meanSignal(counter,k);
        end
    end
end





%%  REMOVE ANY OUTLIERS, ACCORDING TO FLAGS

Hyb_data.medianSignal = Hyb_rawdata.medianSignal;
Hyb_data.meanSignal = Hyb_rawdata.meanSignal;
Assay_data.medianSignal = Assay_rawdata.medianSignal;
Assay_data.meanSignal = Assay_rawdata.meanSignal;
Assay_data.medianBackground = Assay_rawdata.medianBackground;
Assay_data.meanBackground = Assay_rawdata.meanBackground;

Hyb_data.medianSignal_rmOutliers = Hyb_data.medianSignal;
Hyb_data.meanSignal_rmOutliers = Hyb_data.meanSignal;
Assay_data.medianSignal_rmOutliers = Assay_data.medianSignal;
Assay_data.meanSignal_rmOutliers = Assay_data.meanSignal;
Assay_data.medianBackground_rmOutliers = Assay_data.medianBackground;
Assay_data.meanBackground_rmOutliers = Assay_data.meanBackground;



%% remove datapoints based on background signal:
if FLAG_rmOutliersBKGND == 1
    
    for j = 1:size(Hyb_rawdata.medianSignal,2)
        
        Temp = (1:size(Hyb_data.medianSignal_rmOutliers(:,j),1));
        TempH = Temp'./Temp';        
        TempH(...
            Hyb_rawdata.meanBackground(:,j) > ...
            mean(Hyb_rawdata.meanBackground(:,j)) + ...
            FLAG_rmOutliersBKGND_sigma * ...
            std(Hyb_rawdata.meanBackground(:,j)) ) = NaN;

        Temp2 = (1:size(Assay_data.medianSignal_rmOutliers(:,j),1));
        TempA = Temp2'./Temp2';        
        TempA(...
            Assay_rawdata.meanBackground(:,j) > ...
            mean(Assay_rawdata.meanBackground(:,j)) + ...
            FLAG_rmOutliersBKGND_sigma * ...
            std(Assay_rawdata.meanBackground(:,j)) ) = NaN;
        
        Hyb_data.medianSignal_rmOutliers(:,j) = ...
            Hyb_data.medianSignal_rmOutliers(:,j).*TempH;
        Hyb_data.meanSignal_rmOutliers(:,j) = ...
            Hyb_data.meanSignal_rmOutliers(:,j).*TempH;
        
        Assay_data.medianSignal_rmOutliers(:,j) = ...
            Assay_data.medianSignal_rmOutliers(:,j).*TempA;
        Assay_data.meanSignal_rmOutliers(:,j) = ...
            Assay_data.meanSignal_rmOutliers(:,j).*TempA;
        Assay_data.medianBackground_rmOutliers(:,j) = ...
            Assay_data.medianBackground_rmOutliers(:,j).*TempA;
        Assay_data.meanBackground_rmOutliers(:,j) = ...
            Assay_data.meanBackground_rmOutliers(:,j).*TempA;
        
    end
    clear Temp TempH Temp2 TempA
end
% Count how many outliers were removed per subarray for Hyb and Assay
Hyb_data.N_rmOutliersBKGND = sum(isnan(Hyb_data.medianSignal_rmOutliers));
Assay_data.N_rmOutliersBKGND = sum(isnan(Assay_data.medianSignal_rmOutliers));






%% Remove based on spot morphology:
if FLAG_rmOutliersMorphology == 1
    
    for j = 1:size(Hyb_rawdata.medianSignal,2)
        
        Temp = (1:size(Hyb_data.medianSignal_rmOutliers(:,j),1));
        TempH = Temp'./Temp';        
        TempH(...
            ( Hyb_rawdata.stdevSignal(:,j)./Hyb_rawdata.meanSignal(:,j)) > ...
            mean( Hyb_rawdata.stdevSignal(:,j)./Hyb_rawdata.meanSignal(:,j) ) + ...
            FLAG_rmOutliersMorphology_sigma * ...
            std( Hyb_rawdata.stdevSignal(:,j)./Hyb_rawdata.meanSignal(:,j) ) ...
            ) = NaN;

        Temp2 = (1:size(Hyb_data.medianSignal_rmOutliers(:,j),1));
        TempA = Temp2'./Temp2';        
        TempA(...
            ( Assay_rawdata.stdevSignal(:,j)./Assay_rawdata.meanSignal(:,j)) > ...
            mean( Assay_rawdata.stdevSignal(:,j)./Assay_rawdata.meanSignal(:,j) ) + ...
            FLAG_rmOutliersMorphology_sigma * ...
            std( Assay_rawdata.stdevSignal(:,j)./Assay_rawdata.meanSignal(:,j) ) ...
            ) = NaN;

        
        Hyb_data.medianSignal_rmOutliers(:,j) = ...
            Hyb_data.medianSignal_rmOutliers(:,j).*TempH;
        Hyb_data.meanSignal_rmOutliers(:,j) = ...
            Hyb_data.meanSignal_rmOutliers(:,j).*TempH;
        
        Assay_data.medianSignal_rmOutliers(:,j) = ...
            Assay_data.medianSignal_rmOutliers(:,j).*TempA;
        Assay_data.meanSignal_rmOutliers(:,j) = ...
            Assay_data.meanSignal_rmOutliers(:,j).*TempA;
        Assay_data.medianBackground_rmOutliers(:,j) = ...
            Assay_data.medianBackground_rmOutliers(:,j).*TempA;
        Assay_data.meanBackground_rmOutliers(:,j) = ...
            Assay_data.meanBackground_rmOutliers(:,j).*TempA;
        
    end
    clear Temp TempH Temp2 TempA
end
% Count how many outliers were removed per subarray for Hyb and Assay
Hyb_data.N_rmOutliersMorphology = sum(isnan(Hyb_data.medianSignal_rmOutliers));
Assay_data.N_rmOutliersMorphology = sum(isnan(Assay_data.medianSignal_rmOutliers));








%% SORT MICROARRAY DATA
%   use left and right indices from the .gal file:
%   Work with rmOutlier datasets from now on


% Sort for median values
tempHyb = Hyb_data.medianSignal_rmOutliers;
tempAssay = Assay_data.medianSignal_rmOutliers;
Hyb_rmOutliers_Sorted = [ ];
Assay_rmOutliers_Sorted = [ ];
for i = [1 2 3 4 5 6]
    if i == 2 | i == 4 | i == 6
        temp2Hyb = tempHyb(sort_indices_Right,i);
        temp2Assay = tempAssay(sort_indices_Right,i);
    else
        temp2Hyb = tempHyb(sort_indices_Left,i);
        temp2Assay = tempAssay(sort_indices_Left,i);
    end
    Hyb_rmOutliers_Sorted = [Hyb_rmOutliers_Sorted temp2Hyb];
    Assay_rmOutliers_Sorted = [ Assay_rmOutliers_Sorted temp2Assay];
end
Hyb_data.medianSignal_rmOutliers_Sorted = Hyb_rmOutliers_Sorted;
Assay_data.medianSignal_rmOutliers_Sorted = Assay_rmOutliers_Sorted;
clear tempHyb tempAssay


% Sort for mean values
tempHyb = Hyb_data.meanSignal_rmOutliers;
tempAssay = Assay_data.meanSignal_rmOutliers;
Hyb_rmOutliers_Sorted = [ ];
Assay_rmOutliers_Sorted = [ ];
for i = [1 2 3 4 5 6]
    if i == 2 | i == 4 | i == 6
        temp2Hyb = tempHyb(sort_indices_Right,i);
        temp2Assay = tempAssay(sort_indices_Right,i);
    else
        temp2Hyb = tempHyb(sort_indices_Left,i);
        temp2Assay = tempAssay(sort_indices_Left,i);
    end
    Hyb_rmOutliers_Sorted = [Hyb_rmOutliers_Sorted temp2Hyb];
    Assay_rmOutliers_Sorted = [ Assay_rmOutliers_Sorted temp2Assay];
end
Hyb_data.meanSignal_rmOutliers_Sorted = Hyb_rmOutliers_Sorted;
Assay_data.meanSignal_rmOutliers_Sorted = Assay_rmOutliers_Sorted;
clear tempHyb tempAssay



% Sort for assay background values
tempAssay1 = Assay_data.medianBackground_rmOutliers;
tempAssay2 = Assay_data.meanBackground_rmOutliers;
Assay_rmOutliers_Sorted1 = [ ];
Assay_rmOutliers_Sorted2 = [ ];

for i = [1 2 3 4 5 6]
    if i == 2 | i == 4 | i == 6
        temp2Assay1 = tempAssay1(sort_indices_Right,i);
        temp2Assay2 = tempAssay2(sort_indices_Right,i);
    else
        temp2Assay1 = tempAssay1(sort_indices_Left,i);
        temp2Assay2 = tempAssay2(sort_indices_Left,i);

    end
    Assay_rmOutliers_Sorted1 = [ Assay_rmOutliers_Sorted1 temp2Assay1];
    Assay_rmOutliers_Sorted2 = [ Assay_rmOutliers_Sorted2 temp2Assay2];
end
Assay_data.medianBackground_rmOutliers_Sorted = Assay_rmOutliers_Sorted1;
Assay_data.meanBackground_rmOutliers_Sorted = Assay_rmOutliers_Sorted2;
clear tempAssay1 tempAssay2





%% COMPRESS MICROARRAY DATA INTO SINGLE DATA POINTS FOR EACH ACE

for j=1:N_SubArrays
    
    counter = 0; %reset count within each subarray
    
    for i = 1:N_ReplicatesPerACE:N_ACEs*N_ReplicatesPerACE
        
        counter = counter+1; %count forward by N_ReplicatesPerACE
        
        % average over N_ReplicatesPerACE for each ACE
        Hyb_data.medianSRS_C(counter,j) = nanmean(Hyb_data.medianSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Hyb_data.medianSRS_C_Std(counter,j) = nanstd(Hyb_data.medianSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Hyb_data.meanSRS_C(counter,j) = nanmean(Hyb_data.meanSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Hyb_data.meanSRS_C_Std(counter,j) = nanstd(Hyb_data.meanSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Hyb_data.SRS_C_NaNcount(counter,j) = sum(isnan(Hyb_data.medianSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j)));
        
        Assay_data.medianSRS_C(counter,j) = nanmean(Assay_data.medianSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Assay_data.medianSRS_C_Std(counter,j) = nanstd(Assay_data.medianSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Assay_data.meanSRS_C(counter,j) = nanmean(Assay_data.meanSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Assay_data.meanSRS_C_Std(counter,j) = nanstd(Assay_data.meanSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Assay_data.SRS_C_NaNcount(counter,j) = sum(isnan(Assay_data.medianSignal_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j)));
        
        Assay_data.medianBRS_C(counter,j) = nanmean(Assay_data.medianBackground_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Assay_data.meanBRS_C(counter,j) = nanmean(Assay_data.meanBackground_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Assay_data.medianBRS_C_std(counter,j) = nanstd(Assay_data.medianBackground_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        Assay_data.meanBRS_C_std(counter,j) = nanstd(Assay_data.meanBackground_rmOutliers_Sorted(i:i+N_ReplicatesPerACE-1,j));
        
        
    end
end


%% Remove any ACEs that contain too few replicate conditions:
Hyb_data.medianSRS_C(Hyb_data.SRS_C_NaNcount > ...
    N_ReplicatesPerACE - minACEsAfterOutlierRemoval) = NaN;
Hyb_data.meanSRS_C(Hyb_data.SRS_C_NaNcount > ...
    N_ReplicatesPerACE - minACEsAfterOutlierRemoval) = NaN;


Assay_data.medianSRS_C(Assay_data.SRS_C_NaNcount > ...
    N_ReplicatesPerACE - minACEsAfterOutlierRemoval) = NaN;
Assay_data.meanSRS_C(Assay_data.SRS_C_NaNcount > ...
    N_ReplicatesPerACE - minACEsAfterOutlierRemoval) = NaN;

Assay_data.medianBRS_C(Assay_data.SRS_C_NaNcount > ...
    N_ReplicatesPerACE - minACEsAfterOutlierRemoval) = NaN;
Assay_data.meanBRS_C(Assay_data.SRS_C_NaNcount > ...
    N_ReplicatesPerACE - minACEsAfterOutlierRemoval) = NaN;




%% REORDER ARRAY DATA FROM HIGH TO LOW TARGET CONCENTRATION
% Reorder blocks from high to low target, Buffer, Blank conditions


Hyb_data.medianSRS_C_Reorder = Hyb_data.medianSRS_C(:,reorderBlocks);
Hyb_data.meanSRS_C_Reorder = Hyb_data.meanSRS_C(:,reorderBlocks);

Assay_data.medianSRS_C_Reorder = Assay_data.medianSRS_C(:,reorderBlocks);
Assay_data.meanSRS_C_Reorder = Assay_data.meanSRS_C(:,reorderBlocks);

Assay_data.medianBRS_C_Reorder = Assay_data.medianBRS_C(:,reorderBlocks);
Assay_data.meanBRS_C_Reorder = Assay_data.meanBRS_C(:,reorderBlocks);




%% CALCULATE ACE SCORES:

ACE_Scores.median = Assay_data.medianSRS_C(:,reorderBlocks)./Hyb_data.medianSRS_C(:,reorderBlocks);
ACE_Scores.median_Std = ...
    ACE_Scores.median .* sqrt( ...
    ( Hyb_data.medianSRS_C_Std(:,reorderBlocks)./Hyb_data.medianSRS_C(:,reorderBlocks) ).^2 + ... 
    ( Assay_data.medianSRS_C_Std(:,reorderBlocks)./Assay_data.medianSRS_C(:,reorderBlocks) ).^2 ...
    );


ACE_Scores.mean = Assay_data.meanSRS_C(:,reorderBlocks)./Hyb_data.meanSRS_C(:,reorderBlocks);
ACE_Scores.mean_Std = ...
    ACE_Scores.mean .* sqrt( ...
    (Hyb_data.meanSRS_C_Std(:,reorderBlocks)./Hyb_data.meanSRS_C(:,reorderBlocks)).^2 + ... 
    (Assay_data.meanSRS_C_Std(:,reorderBlocks)./Assay_data.meanSRS_C(:,reorderBlocks)).^2 ...
    );

% ACE_Scores.backgroundMedian = Assay_data.medianBRS_C(:,reorderBlocks);
% ACE_Scores.backgroundMean = Assay_data.meanBRS_C(:,reorderBlocks);




%% NORMALIZE ARRAY DATA TO BLANK CONDITION (BLANK = 1 for every ACE)

for j = 1:size(ACE_Scores.median,2)
    
    ACE_Scores.median_Normalized(:,j) = ACE_Scores.median(:,j)./ ACE_Scores.median(:,end);
    
    ACE_Scores.median_Normalized_std(:,j) = ...
        ACE_Scores.median_Normalized(:,j) .* sqrt( ...
        ( ACE_Scores.median_Std(:,j) ./ ACE_Scores.median(:,j) ).^2 + ...
        ( ACE_Scores.median_Std(:,end) ./ ACE_Scores.median(:,end) ).^2 ...
        );
    
	ACE_Scores.mean_Normalized(:,j) = ACE_Scores.mean(:,j)./ ACE_Scores.mean(:,end);
    
    ACE_Scores.mean_Normalized_std(:,j) = ...
        ACE_Scores.mean_Normalized(:,j) .* sqrt( ...
        ( ACE_Scores.mean_Std(:,j) ./ ACE_Scores.mean(:,j) ).^2 + ...
        ( ACE_Scores.mean_Std(:,end) ./ ACE_Scores.mean(:,end) ).^2 ...
        );
    
end





%% REMOVE POORLY HYBRIDIZED ACEs FROM RESULTS (based on HYB and/or ASSAY microarray)

ACE_Scores.median_Normalized_rmLowHyb = ACE_Scores.median_Normalized;
ACE_Scores.mean_Normalized_rmLowHyb = ACE_Scores.mean_Normalized;

if FLAG_rmLowHybdata == 1 | FLAG_rmLowHybdata == 2
    ACE_Scores.median_Normalized_rmLowHyb( ...
        Hyb_data.medianSRS_C <= ...
        FLAG_rmLowHybdata_HybFloor) = NaN; 
    
    ACE_Scores.mean_Normalized_rmLowHyb( ...
        Hyb_data.meanSRS_C <= ...
        FLAG_rmLowHybdata_HybFloor) = NaN;
end

if FLAG_rmLowHybdata == 2
    ACE_Scores.median_Normalized_rmLowHyb( ...
        Assay_data.medianSRS_C <= ...
        FLAG_rmLowHybdata_AssayFloor) = NaN; 
    
    ACE_Scores.mean_Normalized_rmLowHyb( ...
        Assay_data.meanSRS_C <= ...
        FLAG_rmLowHybdata_AssayFloor) = NaN; 
end





%% Remove any ACEs that contain an NaN value on any subarray
for i = 1:size(ACE_Scores.median_Normalized_rmLowHyb,1)
    if sum(isnan(ACE_Scores.median_Normalized_rmLowHyb(i,:))) > 0
        ACE_Scores.median_Normalized_rmLowHyb(i,:) = NaN;
    end
    
    if sum(isnan(ACE_Scores.mean_Normalized_rmLowHyb(i,:))) > 0
        ACE_Scores.mean_Normalized_rmLowHyb(i,:) = NaN;
    end    
    
end




%% CALCUALTE MICHAELIS-MENTEN VALUES (Km, Vmax)
% Skip if Flag_SkipKm = 1

if Flag_SkipKm == 1
    
else
warning off
Vmax_min = -0.1;
Vmax_max = 1.2;
[   ACE_Scores.Km, ACE_Scores.Km_Std, ...
    ACE_Scores.Vmax, ACE_Scores.Vmax_Std, ...
    ACE_Scores.MSError] = ...
    extractMichaelisMenten(ACE_Scores.median_Normalized_rmLowHyb(:,1:end-1), LigandConc, Vmax_min, Vmax_max);

warning on
end




%% LOAD ACE SEQUENCES (5' -> 3')

cd('RawDataFiles/');
cd(workingfolder)

warning off bioinfo:oligoprop:SeqLengthTooShort; % Turn off warnings from matlab's oligoprop function

[ACEs.sequences, ACEs.Tm_NN_Matlab, ACEs.deltaG, ACEs.deltaG_self] = ...
    extractOligoProperties(filename_ACEsequences, N_ACEs);

cd ../..




%% CALCUALTE VARIABLES FOR HEATMAP PLOTS

lengthCounterLeft = 1;
lengthCounterRight = 0;
designCounter = 1;
misMatchCounter = 0;
ACEs.plotMatrix = {};


% scan over all ACE sequences defined in the file
for i = 1:size(ACEs.sequences,1)-1
    
    lengthCounterRight = lengthCounterRight+1;
    
    % check if we are in the first mismatch regime
    if i >= misMatchIndex(1) && i <= misMatchIndex(2)+1
        
        misMatchCounter = misMatchCounter+1;
        
        if misMatchCounter == misMatch(1)
            %End of one mismatched ACE 5' location, go to next
            
            ACEs.plotMatrix{designCounter,1} = lengthCounterLeft:lengthCounterRight;
            designCounter = designCounter+1;
            lengthCounterLeft = lengthCounterRight+1;
            misMatchCounter = 0;
        end

    % if not, check if we are in the second mismatch regime
    elseif i >= misMatchIndex(3) && i <= misMatchIndex(4)+1
        
        misMatchCounter = misMatchCounter+1;
        
        if misMatchCounter == misMatch(2)
            % End of one mismatched ACE 5' location, go to next
            
            ACEs.plotMatrix{designCounter,1} = lengthCounterLeft:lengthCounterRight;
            designCounter = designCounter+1;
            lengthCounterLeft = lengthCounterRight+1;
            misMatchCounter = 0;
        end

    elseif  size(ACEs.sequences{i},2) == size(ACEs.sequences{i+1},2) ...
            &&  strncmpi(ACEs.sequences(i+1),'TTTTT',5) == 0        
        
    elseif  strncmpi(ACEs.sequences(i),'TTTTT',5) == 1 ...
            &&  strncmpi(ACEs.sequences(i+1),'TTTTT',5) == 1
        
        % If working with 5' PolyT ACEs, test for size increase
        if  size(ACEs.sequences{i},2) +1 == size(ACEs.sequences{i+1},2)
        
            ACEs.plotMatrix{designCounter,1} = lengthCounterLeft:lengthCounterRight;
            designCounter = designCounter+1;
            lengthCounterLeft = lengthCounterRight+1;
            
        end
        
    elseif  size(ACEs.sequences{i},2) ~= size(ACEs.sequences{i+1},2) ...
            || strncmpi(ACEs.sequences(i+1),'TTTTT',5) == 1
        
        % If size changes, start new matrix of indexes
        
        ACEs.plotMatrix{designCounter,1} = lengthCounterLeft:lengthCounterRight;
        designCounter = designCounter+1;
        lengthCounterLeft = lengthCounterRight+1;
        
        
    end
    
    % Case of last elements in the array, group together:
    if i == size(ACEs.sequences,1)-1
        lengthCounterRight = lengthCounterRight+1;
        ACEs.plotMatrix{designCounter,1} = lengthCounterLeft:lengthCounterRight;
    end
    
end






%% DEFINE ADDITIONAL SCORES FOR MAKING HEAT MAP FIGURES (koff, k*off, ratios)


%remove ACEs that have negative switching...
ACE_Scores.median_Normalized_std(Hyb_data.SRS_C_NaNcount > ...
    N_ReplicatesPerACE - minACEsAfterOutlierRemoval) = NaN;

%Baseline score = Blank - Buffer
ACE_Scores.baselineSwitch = 1-ACE_Scores.median_Normalized_rmLowHyb(:,5);
ACE_Scores.baselineSwitch_std = sqrt(ACE_Scores.median_Normalized_std(:,6).^2 + ACE_Scores.median_Normalized_std(:,5).^2);
ACE_Scores.baselineSwitch_std(isnan(ACE_Scores.baselineSwitch))= NaN;

%MaxSwitch score =  Buffer - Max Ligand
ACE_Scores.maxSwitch = ACE_Scores.median_Normalized_rmLowHyb(:,5)-ACE_Scores.median_Normalized_rmLowHyb(:,1);
ACE_Scores.maxSwitch_std = sqrt(ACE_Scores.median_Normalized_std(:,5).^2 + ACE_Scores.median_Normalized_std(:,1).^2) ;
ACE_Scores.maxSwitch_std(isnan(ACE_Scores.maxSwitch))= NaN;

%SwitchRatio1 score = abs(MaxSwitch)/abs(Baseline)
ACE_Scores.switchRatio1 = abs(ACE_Scores.maxSwitch) ./ abs(ACE_Scores.baselineSwitch) ;
ACE_Scores.switchRatio1_std = ACE_Scores.switchRatio1.*sqrt( (ACE_Scores.maxSwitch_std./ACE_Scores.maxSwitch).^2 + (ACE_Scores.baselineSwitch_std./ACE_Scores.baselineSwitch).^2);

%SwitchRatio2 score = abs[ (MaxSwitch) + (Baseline) ] ./ abs(Baseline)
ACE_Scores.switchRatio2 = abs(ACE_Scores.maxSwitch + ACE_Scores.baselineSwitch)./ abs(ACE_Scores.baselineSwitch) ;
temp = abs(ACE_Scores.maxSwitch + ACE_Scores.baselineSwitch);
temp_std = sqrt(ACE_Scores.maxSwitch_std.^2 + ACE_Scores.baselineSwitch_std.^2);
ACE_Scores.switchRatio2_std = abs(ACE_Scores.switchRatio2).*sqrt((temp_std./temp).^2 + (ACE_Scores.baselineSwitch_std./ACE_Scores.baselineSwitch).^2);
clear temp temp_std




%% LEAST SQUARES FITTING
% Linearly fit rates vs log10 of ligand concentration
for n = 1:size(ACE_Scores.median_Normalized_rmLowHyb,1)
    x = log10(LigandConc);
    p = polyfit(x, ACE_Scores.median_Normalized_rmLowHyb(n,1:end-2),1); % Fit to a line
    slope(n) = p(1);
    intercept(n) = p(2);
    ycalc =  slope(n)*x + intercept(n);
    
    SSres(n) = sum((ACE_Scores.median_Normalized_rmLowHyb(n,1:end-2) - ycalc).^2);
    SStot(n) = sum((ACE_Scores.median_Normalized_rmLowHyb(n,1:end-2) - mean(ACE_Scores.median_Normalized_rmLowHyb(n,1:end-2))).^2);
    
    Rsq(n) = 1 - SSres(n)/SStot(n);
end
ACE_Scores.LSFit_Slope = slope';
ACE_Scores.LSFit_Intercept = intercept';
ACE_Scores.LSFit_Rsq = Rsq';




%% TILE THE DATA INTO HEATMAPS

[ACE_Scores.baselineSwitchTiled5Prime, ACE_Scores.baselineSwitchTiled3Prime] = tileSmoothHeatMapData(...
    ACE_Scores.baselineSwitch, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

[ACE_Scores.baselineSwitch_stdTiled5Prime, ACE_Scores.baselineSwitch_stdTiled3Prime] = tileSmoothHeatMapData(...
    ACE_Scores.baselineSwitch_std, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);


[ACE_Scores.maxSwitchTiled5Prime, ACE_Scores.maxSwitchTiled3Prime] = tileSmoothHeatMapData(...
    ACE_Scores.maxSwitch, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

[ACE_Scores.maxSwitch_stdTiled5Prime, ACE_Scores.maxSwitch_stdTiled3Prime] = tileSmoothHeatMapData(...
    ACE_Scores.maxSwitch_std, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);


[ACE_Scores.switchRatio1Tiled5Prime, ACE_Scores.switchRatio1Tiled3Prime] = tileSmoothHeatMapData(...
    ACE_Scores.switchRatio1, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

[ACE_Scores.switchRatio2Tiled5Prime, ACE_Scores.switchRatio2Tiled3Prime] = tileSmoothHeatMapData(...
    ACE_Scores.switchRatio2, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);


% Also Tile Km and Vmax values
if Flag_SkipKm == 1
else
[ACE_Scores.VmaxTiled5Prime, ACE_Scores.VmaxTiled3Prime] = tileSmoothHeatMapData( ...
    ACE_Scores.Vmax, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

[ACE_Scores.KmTiled5Prime, ACE_Scores.KmTiled3Prime] = tileSmoothHeatMapData(...
    ACE_Scores.Km, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);
end




%% TILE HYB, ASSAY AND INTERARRAY DATA INTO HEATMAPS

% calculare interarray loss of fluorescence in the BLANK condition:
ACE_Scores.interarray = 100 - Gain_Hyb/Gain_Assay*100*Assay_data.medianSRS_C(:,6)./(Hyb_data.medianSRS_C(:,6));
ACE_Scores.interarray(isnan(ACE_Scores.baselineSwitch))= NaN;


[ACE_Scores.interarrayTiled5Prime, ACE_Scores.interarrayTiled3Prime] = tileSmoothHeatMapData(...
    ACE_Scores.interarray, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);


[Hyb_data.medianSRS_C_Tiled5Prime, Hyb_data.medianSRS_C_Tiled3Prime] = tileSmoothHeatMapData(...
    Hyb_data.medianSRS_C(:,1), ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);


[Assay_data.medianSRS_C_Tiled5Prime, Assay_data.medianSRS_C_Tiled3Prime] = tileSmoothHeatMapData(...
    Assay_data.medianSRS_C(:,1), ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);


[ACEs.Tm_NN_Matlab_Tiled5Prime, ACEs.Tm_NN_Matlab_Tiled3Prime] = tileSmoothHeatMapData(...
    ACEs.Tm_NN_Matlab, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

[ACEs.deltaG_Tiled5Prime, ACEs.deltaG_Tiled3Prime] = tileSmoothHeatMapData(...
    ACEs.deltaG, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

[ACEs.deltaG_self_Tiled5Prime, ACEs.deltaG_self_Tiled3Prime] = tileSmoothHeatMapData(...
    ACEs.deltaG_self, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);


% Tile with background subtraction:
[Assay_data.medianBackgroundTiled5Prime, Assay_data.medianBackgroundTiled3Prime] = tileSmoothHeatMapData(...
    Assay_data.medianSRS_C_Reorder(:,1) - Assay_data.medianBRS_C_Reorder(:,1), ...
    ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

[Assay_data.meanBackgroundTiled5Prime, Assay_data.meanBackgroundTiled3Prime] = tileSmoothHeatMapData(...
    Assay_data.meanSRS_C_Reorder(:,1) - Assay_data.meanBRS_C_Reorder(:,1), ...
    ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);




%% In this section, decompose the Heatmaps generated above to only show e.g. 7-15mers, only MM's, etc...

HeatMapDecompCounter = [];
HeatMapDecompCounter(1) = 1;
for i = 1:size(ACEs.plotMatrix, 1)
    for j = ACEs.plotMatrix{i}(1):ACEs.plotMatrix{i}(end)
        if j == misMatchIndex(1)
            HeatMapDecompCounter(2) = i-1;
            HeatMapDecompCounter(3) = i;
        elseif j == misMatchIndex(2)
            HeatMapDecompCounter(4) = i;
        elseif j == misMatchIndex(3)
            HeatMapDecompCounter(5) = i;
        elseif j == misMatchIndex(4)
            HeatMapDecompCounter(6) = i;
        end
    end
end
HeatMapDecompCounter(7) = size(ACEs.plotMatrix, 1);




%% Tile Interarray Data, but separated into perfect match or MisMatch sections

[ACE_Scores_Perf.baselineSwitchTiled5Prime, ACE_Scores_Perf.baselineSwitchTiled3Prime, ACE_Scores_Perf.baselineSwitchTiled5Prime_Smooth, ACE_Scores_Perf.baselineSwitch_Counts] = tileSmoothHeatMapData(...
    ACE_Scores.baselineSwitch, ACEs.sequences, ACEs.plotMatrix(HeatMapDecompCounter(1):HeatMapDecompCounter(2)), startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

[ACE_Scores_Perf.maxSwitchTiled5Prime, ACE_Scores_Perf.maxSwitchTiled3Prime, ACE_Scores_Perf.maxSwitchTiled5Prime_Smooth, ACE_Scores_Perf.maxSwitch_Counts] = tileSmoothHeatMapData(...
    ACE_Scores.maxSwitch, ACEs.sequences, ACEs.plotMatrix(HeatMapDecompCounter(1):HeatMapDecompCounter(2)), startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);


[ACE_Scores_MM1.baselineSwitchTiled5Prime, ACE_Scores_MM1.baselineSwitchTiled3Prime, ACE_Scores_MM1.baselineSwitchTiled5Prime_Smooth, ACE_Scores_MM1.baselineSwitch_Counts] = tileSmoothHeatMapData(...
    ACE_Scores.baselineSwitch, ACEs.sequences, ACEs.plotMatrix(HeatMapDecompCounter(3):HeatMapDecompCounter(4)), startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

[ACE_Scores_MM1.maxSwitchTiled5Prime, ACE_Scores_MM1.maxSwitchTiled3Prime, ACE_Scores_MM1.maxSwitchTiled5Prime_Smooth, ACE_Scores_MM1.maxSwitch_Counts] = tileSmoothHeatMapData(...
    ACE_Scores.maxSwitch, ACEs.sequences, ACEs.plotMatrix(HeatMapDecompCounter(3):HeatMapDecompCounter(4)), startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);



%% Project mismatches onto aptamer sequence
ACE_Scores_MM1.maxSwitch_ProjectOnAptamer = nanmean(ACE_Scores_MM1.maxSwitchTiled5Prime');





%% Quantitative data analysis
if QuantitativeFlag == 1;
    cd('RawDataFiles/');
    cd(workingfolder)
    
    
    % Read in file
    StartRow = 1;
    StartCol = 0;
    % Open files, read file contents
    try
        QuantExtractedData.raw = csvread(filename_ExtractedData,StartRow,StartCol);
    catch me
        disp('Can not find ExtractedData file')
    end
    
    % Remove any ACEs that contain an NaN value on any subarray
    for i = 1:size(QuantExtractedData.raw,1)
        if sum(isnan(QuantExtractedData.raw(i,:))) > 0
            QuantExtractedData.raw(i,:) = NaN;
        end
    end
    
    QuantExtractedData.processed(:,1:8) = QuantExtractedData.raw(:,1:8);
    QuantExtractedData.processed(:,9) = mean(QuantExtractedData.raw(:,8:9),2);
    
    
    % CALCUALTE MICHAELIS-MENTEN VALUES (Km, Vmax)
    
    warning off
    QuantLigandConc = [.01 0.002 .0004 .00008 .000016 .0000032 .00000064 .000000128]; % ATP concentrations, Molar
    
    QuantVmax_min = -0.1;
    QuantVmax_max = 1.2;
    [   QuantExtractedData.Km, QuantExtractedData.Km_Std, ...
        QuantExtractedData.Vmax, QuantExtractedData.Vmax_Std, ...
        QuantExtractedData.MSError] = ...
        extractMichaelisMenten(QuantExtractedData.processed, QuantLigandConc, QuantVmax_min, QuantVmax_max);
    
    warning on
    
    cd ../..
    
end





%% Surface density data analysis for ATP DNA aptamer
if SurfaceDensityFlag == 1;
    
    cd('RawDataFiles/');
    cd(workingfolder)

    % Read in file
    StartRow = 2;
    StartCol = 0;
    % Open files, read file contents
    try
        SurfaceDensity.raw = csvread(filename_SurfaceDensityData,StartRow,StartCol);
    catch me
        disp('Can not find SurfaceDensityData file')
    end

    
    % Remove any ACEs that contain an NaN value on any subarray
    for i = 1:size(SurfaceDensity.raw,1)
        if sum(isnan(SurfaceDensity.raw(i,:))) > 0
            SurfaceDensity.raw(i,:) = NaN;
        end
    end
    
    
    SurfaceDensity.processed(:,1) = SurfaceDensity.raw(:,1);
    SurfaceDensity.processed(:,2) = SurfaceDensity.raw(:,2);
    SurfaceDensity.processed(:,3:8) = SurfaceDensity.raw(:,5:10);
    
    
    SurfaceDensity.baseline(:,1) = 1-SurfaceDensity.processed(:,1);
    SurfaceDensity.baseline(:,2) = 1-SurfaceDensity.processed(:,3);
    SurfaceDensity.baseline(:,3) = 1-SurfaceDensity.processed(:,5);
    SurfaceDensity.baseline(:,4) = 1-SurfaceDensity.processed(:,7);

    SurfaceDensity.maxligand(:,1) = SurfaceDensity.processed(:,1)-SurfaceDensity.processed(:,2);
    SurfaceDensity.maxligand(:,2) = SurfaceDensity.processed(:,3)-SurfaceDensity.processed(:,4);
    SurfaceDensity.maxligand(:,3) = SurfaceDensity.processed(:,5)-SurfaceDensity.processed(:,6);
    SurfaceDensity.maxligand(:,4) = SurfaceDensity.processed(:,7)-SurfaceDensity.processed(:,8);
   
    cd ../..
end





%% TBA RED and GREEN channel Analysis
if CompareRedGreenChannelsFlag == 1;
    
    cd('RawDataFiles/');
    cd(workingfolder)
    
    
    % Read in file
    StartRow = 1;
    StartCol = 0;
    % Open files, read file contents
    try
        data_REDGREEN.raw = csvread(filename_RED_GREEN_data,StartRow,StartCol);
    catch me
        disp('Can not find RED_GREEN_Data file')
    end

    
    data_REDGREEN.HybGREEN = data_REDGREEN.raw(:,1);
    data_REDGREEN.Assay1GREEN = data_REDGREEN.raw(:,2);
    data_REDGREEN.Assay2GREEN = data_REDGREEN.raw(:,3);
    data_REDGREEN.Assay1RED = data_REDGREEN.raw(:,4);
    data_REDGREEN.Assay2RED = data_REDGREEN.raw(:,5);
    
    data_REDGREEN.NormHyb1 = data_REDGREEN.Assay1RED ./ data_REDGREEN.HybGREEN;
    data_REDGREEN.NormHyb2 = data_REDGREEN.Assay2RED ./ data_REDGREEN.HybGREEN;
    data_REDGREEN.NormAssay1 = data_REDGREEN.Assay1RED ./ data_REDGREEN.Assay1GREEN;
    data_REDGREEN.NormAssay2 = data_REDGREEN.Assay2RED ./ data_REDGREEN.Assay2GREEN;

    [data_REDGREEN.AbsoluteAssay1Tiled5Prime, ACE_Scores_Perf.NormAssay1Tiled3Prime, ~, ~ ] = tileSmoothHeatMapData(...
    data_REDGREEN.Assay1RED, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);
    
    [data_REDGREEN.NormAssay1Tiled5Prime, ACE_Scores_Perf.NormAssay1Tiled3Prime, ~, ~ ] = tileSmoothHeatMapData(...
    data_REDGREEN.NormAssay1, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

    [data_REDGREEN.NormAssay2Tiled5Prime, ACE_Scores_Perf.NormAssay2Tiled3Prime, ~, ~ ] = tileSmoothHeatMapData(...
    data_REDGREEN.NormAssay2, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

    cd ../..

end




%% TBA Sodium vs. Potassium Analysis
if CompareSodiumPotassiumFlag == 1;

    cd('RawDataFiles/');
    cd(workingfolder)
    
    
    % Read in file
    StartRow = 1;
    StartCol = 0;
    % Open files, read file contents
    try
        data_SodiumPotassium.raw = csvread(filename_Sodium_Potassium_data,StartRow,StartCol);
    catch me
        disp('Can not find Sodium vs Potassium data file')
    end

    
    data_SodiumPotassium.Sodium.Baseline = data_SodiumPotassium.raw(:,1);
    data_SodiumPotassium.Sodium.MaxLigand = data_SodiumPotassium.raw(:,2);
    data_SodiumPotassium.Potassium.Baseline = data_SodiumPotassium.raw(:,3);
    data_SodiumPotassium.Potassium.MaxLigand = data_SodiumPotassium.raw(:,4);

    cd ../..
    
end





%% Riboswitch covariation data analysis
if RiboswitchFlagCovariation == 1
    
    cd('RawDataFiles/');
    cd(workingfolder)
    
    
    % Read in file
    StartRow = 1;
    StartCol = 2;
    % Open files, read file contents
    try
        Covariation.raw = csvread(filename_CovationData,StartRow,StartCol);
    catch me
        disp('Can not find CovariationData file')
    end
    
    
    Covariation.P1 = [ Covariation.raw(6:12,:) ; Covariation.raw(66:72,:)];
    Covariation.J12 = Covariation.raw(13:15,:);
    Covariation.P2 = [ Covariation.raw(16:22,:) ; Covariation.raw(30:36,:)];
    Covariation.L2 = Covariation.raw(23:29,:);
    Covariation.J23 = Covariation.raw(37:44,:);
    Covariation.P3 = [ Covariation.raw(45:50,:) ; Covariation.raw(58:63,:)];
    Covariation.L3 = Covariation.raw(51:57,:);
    Covariation.J31 = Covariation.raw(64:65,:);
    Covariation.P4 = Covariation.raw(72:78,:);
    Covariation.P5 = Covariation.raw(79:87,:);
    
    cd ../..
    
end




%% Average Data for replicate riboswitch concentrations

if RiboswitchFlag == 1
    
% % 22C 4 mM Mg trial
ACE_Scores_Average.data(:,1) = (ACE_Scores.median_Normalized_rmLowHyb(:,1) + ACE_Scores.median_Normalized_rmLowHyb(:,2))/2;
ACE_Scores_Average.data(:,2) = (ACE_Scores.median_Normalized_rmLowHyb(:,4) + ACE_Scores.median_Normalized_rmLowHyb(:,5))/2;
ACE_Scores_Average.data(:,3) = (ACE_Scores.median_Normalized_rmLowHyb(:,6) + ACE_Scores.median_Normalized_rmLowHyb(:,6))/2;

elseif RiboswitchFlag == 2
    
% 10C 12 and 0 mM, and 23C 12 mM
ACE_Scores_Average.data(:,1) = (ACE_Scores.median_Normalized_rmLowHyb(:,1) + ACE_Scores.median_Normalized_rmLowHyb(:,2))/2;
ACE_Scores_Average.data(:,2) = (ACE_Scores.median_Normalized_rmLowHyb(:,3) + ACE_Scores.median_Normalized_rmLowHyb(:,4))/2;
ACE_Scores_Average.data(:,3) = (ACE_Scores.median_Normalized_rmLowHyb(:,5) + ACE_Scores.median_Normalized_rmLowHyb(:,6))/2;

end



if RiboswitchFlag == 1 || RiboswitchFlag == 2

ACE_Scores_Average.Normalized(:,1) = ACE_Scores_Average.data(:,1)./ACE_Scores_Average.data(:,3);
ACE_Scores_Average.Normalized(:,2) = ACE_Scores_Average.data(:,2)./ACE_Scores_Average.data(:,3);
ACE_Scores_Average.Normalized(:,3) = ACE_Scores_Average.data(:,3)./ACE_Scores_Average.data(:,3);

[ACE_Scores_Average.baselineSwitchTiled5Prime] = tileSmoothHeatMapData(...
    1-ACE_Scores_Average.Normalized(:,2), ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);
[ACE_Scores_Average.maxSwitchTiled5Prime] = tileSmoothHeatMapData(...
    ACE_Scores_Average.Normalized(:,2)-ACE_Scores_Average.Normalized(:,1), ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);

ACE_Scores_Average.baselineSwitchTiled5Prime(isnan(ACE_Scores.mean_Normalized_rmLowHyb(1,:)))= NaN;
ACE_Scores_Average.maxSwitchTiled5Prime(isnan(ACE_Scores.mean_Normalized_rmLowHyb(1,:)))= NaN;


figure; sc(ACE_Scores_Average.baselineSwitchTiled5Prime, 'jet', [0, .75], [1,1,1])
figure; sc(ACE_Scores_Average.maxSwitchTiled5Prime, 'jet', [0, .20], [1,1,1])

end




%% SAVE all of the data in the results folder, then come back to the parent


cd('CompiledResults/')
cd(workingfolder)

save(strcat('analyzedMicroarrayData', datestr(now), '.mat'))
cd ../..






%% SAVE all of the heat map figures in the results folder, then come back to the parent

cd('CompiledResults/')
cd(workingfolder)


xxx = figure;
hold on
sc((ACE_Scores.interarrayTiled5Prime), [WashingBlank1 WashingBlank2], 'jet', [1 1 1])
title('Loss from washing Blank [%]')
colorbar
hold off
saveas(xxx,'WashBlankLossHeatMap.png')


xxx = figure;
hold on
sc((Hyb_data.medianSRS_C_Tiled5Prime), [HybLim1 HybLim2], 'jet', [1 1 1])
title('After Hyb: Hybridization of Aptamer to Array [RFU]')
colorbar
hold off
saveas(xxx,'HybDataSignalHeatMap.png')


xxx = figure;
hold on
sc((Assay_data.medianSRS_C_Tiled5Prime), [HybLim1 HybLim2], 'jet', [1 1 1])
title('After Assay: Hybridization of Aptamer to Array [RFU]')
colorbar
hold off
saveas(xxx,'AssayDataSignalHeatMap.png')


xxx = figure;
hold on
sc((ACE_Scores.baselineSwitchTiled5Prime*100), [Baseline1 Baseline2], 'jet', [1 1 1])
title('Loss from Buffer only [%]')
colorbar
hold off
saveas(xxx,'Baseline_BufferLoss.png')


xxx = figure;
hold on
sc((ACE_Scores.baselineSwitch_stdTiled5Prime*100), [Baseline1 Baseline2], 'jet', [1 1 1])
title('StDev of Loss from Buffer only [%]')
colorbar
hold off
saveas(xxx,'Baseline_BufferLoss_Stdev.png')


xxx = figure;
hold on
sc((ACE_Scores.maxSwitchTiled5Prime*100), [Ligand1 Ligand2], 'jet', [1 1 1])
title('Loss from max Ligand concentration [%]')
colorbar
hold off
saveas(xxx,'LigandmaxSwitchRate.png')


xxx = figure;
hold on
sc((ACE_Scores.maxSwitch_stdTiled5Prime*100), [Ligand1 Ligand2], 'jet', [1 1 1])
title('StDev of Loss from max Ligand concentration [%]')
colorbar
hold off
saveas(xxx,'LigandmaxSwitchRate_Stdev.png')


xxx = figure;
hold on
sc((ACE_Scores.switchRatio1Tiled5Prime), [Score1_1 Score1_2], 'jet', [1 1 1])
title('Ratio of Max switch over baseline switch (SBR) [fold change]')
colorbar
hold off
saveas(xxx,'ACE_Score_1.png')


xxx = figure;
hold on
sc((ACE_Scores.switchRatio2Tiled5Prime), [Score2_1 Score2_2], 'jet', [1 1 1])
title('Ratio of (Max switch plus Baseline) over baseline switch (SplusBtoBR) [fold change]')
colorbar
hold off
saveas(xxx,'ACE_Score_2.png')


xxx = figure;
hold on
sc((ACE_Scores_MM1.maxSwitch_ProjectOnAptamer), [0.2 0.3], 'jet', [1 1 1])
title(' Max Ligand Switch rate of Mismatches1 Projected on aptamer sequence')
colorbar
hold off
saveas(xxx,'LigandmaxSwitchRate_Mismatch1_Projected.png')


cd ../..






%% SI Figures

cd('CompiledResults/')
cd(workingfolder)

% Put Ceiling on delta G values of ACE Self-Hyb
ACEs.deltaG_self(ACEs.deltaG_self>1.5) = 1.5;


xxx = figure;
hold on
sc((ACEs.Tm_NN_Matlab_Tiled5Prime), [0 80], 'jet', [1 1 1])
title('Tm (Degrees C)')
colorbar
hold off
saveas(xxx,'TmHeatMap.png')


xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
hold on
colormap 'cool'
scatter(ACEs.deltaG,mean(Hyb_data.medianSRS_C,2), 25, ACEs.deltaG_self, 'filled')
title('compare deltaG vs Hyb Signal')
box on
colorbar
xlim([-55, -5])
caxis([-1.5, 1.5])
hold off
saveas(xxx,'deltaG_vs_KHyb_V1.eps')


xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
hold on
colormap 'cool'
scatter(ACEs.deltaG(352:end),mean(Hyb_data.medianSRS_C(352:end,:),2), 25, [.8 .8 .8] , 'filled')
scatter(ACEs.deltaG(1:351),mean(Hyb_data.medianSRS_C(1:351,:),2), 25, ACEs.deltaG_self(1:351), 'filled')
title('compare deltaG vs Hyb Signal')
box on
colorbar
xlim([-55, -5])
caxis([-1.5, 1.5])
hold off
saveas(xxx,'deltaG_vs_KHyb_V2.eps')


xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
hold on
colormap 'cool'
% scatter(ACEs.deltaG(352:end),ACE_Scores.baselineSwitch(352:end), 25, [.8 .8 .8] , 'filled')
scatter(ACEs.deltaG(1:end),ACE_Scores.baselineSwitch(1:end), 25, ACEs.deltaG_self(1:end), 'filled')
title('compare deltaG vs baseline switch rate')
box on
colorbar
xlim([-55, -5])
ylim([0, 0.7])
caxis([-1.5, 1.5])
hold off
saveas(xxx,'deltaG_vs_BufferOnly.eps')


xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
hold on
colormap 'cool'
% scatter(ACEs.deltaG(352:end),ACE_Scores.maxSwitch(352:end), 25, [.5 .5 .5] , 'filled')
scatter(ACEs.deltaG(1:end),ACE_Scores.maxSwitch(1:end), 25, ACEs.deltaG_self(1:end), 'filled')
title('compare deltaG vs baseline switch rate')
box on
colorbar
xlim([-55, -5])
ylim([0, 0.7])
caxis([-1.5, 1.5])
hold off
saveas(xxx,'deltaG_vs_MaxLigand.eps')


xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
hold on
colormap 'cool'
% scatter(ACEs.deltaG(352:end),ACE_Scores.maxSwitch(352:end), 25, [.5 .5 .5] , 'filled')
scatter(ACE_Scores.baselineSwitch(1:end),ACE_Scores.maxSwitch(1:end), 25, ACEs.deltaG_self(1:end), 'filled')
title('compare baseline switch rate vs max ligand switch rate')
box on
colorbar
xlim([0, 0.7])
ylim([0, 0.7])
caxis([-1.5, 1.5])
hold off
saveas(xxx,'baseline_vs_MaxLigand.eps')




%% Make quantitative figures
if QuantitativeFlag == 1
    
    [QuantExtractedData.Vmax_Tiled5Prime, QuantExtractedData.Vmax_Tiled3Prime] = tileSmoothHeatMapData(...
        QuantExtractedData.Vmax, ACEs.sequences, ACEs.plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex);
    
    QuantExtractedData.Km(isnan(QuantExtractedData.Km)) = 1;
    QuantExtractedData.Km(0.434*QuantExtractedData.Km_Std./QuantExtractedData.Km > 1) = 1;

    
    xxx = figure;
    hold on
    sc((QuantExtractedData.Vmax_Tiled5Prime), [0 0.7], 'jet', [1 1 1])
    title('Vmax (%)')
    colorbar
    hold off
    saveas(xxx,'Vmax_Heatmap.png')
    
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    %9mer
    errorbar(log10(QuantExtractedData.Km(52:75)), 0.434*QuantExtractedData.Km_Std(52:75)./QuantExtractedData.Km(52:75), 'bo', 'MarkerFaceColor', 'b')
    %15mer
%     errorbar(log10(QuantExtractedData.Km(181:198)), 0.434*QuantExtractedData.Km_Std(181:198)./QuantExtractedData.Km(181:198), 'ro', 'MarkerFaceColor', 'r')
    %12mer
    errorbar(log10(QuantExtractedData.Km(121:141)), 0.434*QuantExtractedData.Km_Std(121:141)./QuantExtractedData.Km(121:141), 'gs', 'MarkerFaceColor', 'g')
    title('Quantitative Km Data for 9mer and 15mer ATP DNA DAs')
    box on
    % xlim([0, 0.7])
    % ylim([0, 0.7])
    % caxis([-1.5, 1.5])
    hold off
    saveas(xxx,'ATP_Quantitative_Km.eps')
    
    
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    scatter(1-QuantExtractedData.processed(27:198,9), 1-QuantExtractedData.processed(258:429,9), 25, 'r', 'filled')
    scatter(1-QuantExtractedData.processed(27:198,9), 1-QuantExtractedData.processed(463:634,9), 25, 'b', 'filled')
    title('Compare T10 extension effect, buffer only')
    box on
    legend('5prime T10', '3prime T10')
    xlim([-0.05, 0.3])
    ylim([-0.05, 0.3])
    % caxis([-1.5, 1.5])
    lsline
    pbaspect([1 1 1])
    hold off
    saveas(xxx,'5prime3prime_extensions_Vmin.eps')
    
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    scatter(QuantExtractedData.Vmax(27:198), QuantExtractedData.Vmax(258:429), 25, 'r', 'filled')
    scatter(QuantExtractedData.Vmax(27:198), QuantExtractedData.Vmax(463:634), 25, 'b', 'filled')
    title('Compare T10 extension effect, max ligand')
    box on
    legend('5prime T10', '3prime T10')
    xlim([0, 0.7])
    ylim([0, 0.7])
    % caxis([-1.5, 1.5])
    lsline
    pbaspect([1 1 1])
    hold off
    saveas(xxx,'5prime3prime_extensions_Vmax.eps')
    
end




%% Make surface density figures
if SurfaceDensityFlag == 1;
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    scatter(SurfaceDensity.baseline(1:231,1), SurfaceDensity.baseline(1:231,4), 25, 'm', 'filled')
    scatter(SurfaceDensity.baseline(1:231,1), SurfaceDensity.baseline(1:231,3), 25, 'g', 'filled')
    scatter(SurfaceDensity.baseline(1:231,1), SurfaceDensity.baseline(1:231,2), 25, 'r', 'filled')
    title('Compare surface densities, buffer only')
    box on
    legend('10^9', '10^10', '10^11')

    % caxis([-1.5, 1.5])
    lsline
    pbaspect([1 1 1])
    xlim([-.1, .7])
    ylim([-.1, .7])
    hold off
    saveas(xxx,'SurfaceDensities_baseline.eps')
    
    
        
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    scatter(SurfaceDensity.maxligand(1:231,1), SurfaceDensity.maxligand(1:231,4), 25, 'm', 'filled')
    scatter(SurfaceDensity.maxligand(1:231,1), SurfaceDensity.maxligand(1:231,3), 25, 'g', 'filled')
    scatter(SurfaceDensity.maxligand(1:231,1), SurfaceDensity.maxligand(1:231,2), 25, 'r', 'filled')
    title('Compare surface densities, buffer only')
    box on
    legend('10^9', '10^10', '10^11')
    
    % caxis([-1.5, 1.5])
    lsline
    xlim([-.1, .7])
    ylim([-.1, .7])
    pbaspect([1 1 1])
    hold off
    saveas(xxx,'SurfaceDensities_maxligand.eps')
end




%% Make figures for riboswitches
if RiboswitchFlag == 1 || RiboswitchFlag == 2;
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    scatter(ACEs.deltaG,mean(Hyb_data.medianSRS_C,2), 25, ACEs.deltaG_self, 'filled')
    title('compare deltaG vs Hyb Signal')
    box on
    colorbar
    xlim([-30, 0])
    caxis([-1.5, 1.5])
    pbaspect([1 1 1])
    hold off
    saveas(xxx,'deltaG_vs_KHyb_V1_riboswitch.eps')
    
    
    xxx = figure;
    hold on
    sc(ACE_Scores_Average.baselineSwitchTiled5Prime, [0, .5], 'jet',  [1,1,1])
    title('Buffer only switching of add riboswitch')
    colorbar
    hold off
    saveas(xxx,'Buffer_Only_riboswitch.png')
    
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    scatter(ACEs.deltaG(1:553),mean(Hyb_data.medianSRS_C(1:553,:),2), 25, ACEs.deltaG_self(1:553), 'filled');
    title('compare deltaG vs Hyb Signal')
    box on
    colorbar
    xlim([-25, 0])
    caxis([-1.5, 1.5])
    pbaspect([1 1 1])
    hold off
    saveas(xxx,'deltaG_vs_KHyb_V1_riboswitch.eps')
    
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    scatter(mean(Hyb_data.medianSRS_C(1:553,:),2), ACE_Scores_Average.Normalized(1:553,2) - ACE_Scores_Average.Normalized(1:553,1), 25, ACEs.deltaG_self(1:553), 'filled');
    title('compare Hyb Signal vs max ligand')
    box on
    colorbar
    xlim([0, 30000])
    caxis([-1.5, 1.5])
    pbaspect([1 1 1])
    hold off
    saveas(xxx,'KHyb_vs_maxLigand_riboswitch.eps')
    
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    % scatter(ACEs.deltaG(352:end),ACE_Scores.maxSwitch(352:end), 25, [.5 .5 .5] , 'filled')
    scatter(1-ACE_Scores_Average.Normalized(1:553,2),ACE_Scores_Average.Normalized(1:553,2) - ACE_Scores_Average.Normalized(1:553,1), 25, ACEs.deltaG_self(1:553), 'filled')
    title('compare baseline switch rate vs max ligand switch rate')
    box on
    colorbar
    xlim([0, 0.5])
    ylim([-0.05, 0.25])
    caxis([-1.5, 1.5])
    pbaspect([1 1 1])
    hold off
    saveas(xxx,'baseline_vs_MaxLigand_riboswitch.eps')
    
    
    % plot covariation data
    if RiboswitchFlagCovariation == 1
        
        xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
        hold on
        colormap 'cool'
        scatter(Covariation.P1(:,1),Covariation.P1(:,2), 75, 'bo', 'filled');
        scatter(Covariation.P2(:,1),Covariation.P2(:,2), 75, 'rs', 'filled');
        scatter(Covariation.P3(:,1),Covariation.P3(:,2), 75, 'g^', 'filled');
        scatter(Covariation.P4(:,1),Covariation.P4(:,2), 75, 'md', 'filled');
        scatter(Covariation.P5(:,1),Covariation.P5(:,2), 75, 'cv', 'filled');
        legend('P1', 'P2', 'P3', 'P4', 'P5')
        title('Map 9mers onto riboswitch')
        box on
        % colorbar
        % xlim([-30, 0])
        % caxis([-1.5, 1.5])
        pbaspect([1 1 1])
        hold off
        saveas(xxx,'covariation_9mers_riboswitch_duplexes.eps')
        
        
        xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
        hold on
        colormap 'cool'
        scatter(Covariation.J12(:,1),Covariation.J12(:,2), 75, 'bo', 'filled');
        scatter(Covariation.L2(:,1),Covariation.L2(:,2), 75, 'rs', 'filled');
        scatter(Covariation.J23(:,1),Covariation.J23(:,2), 75, 'g^', 'filled');
        scatter(Covariation.L3(:,1),Covariation.L3(:,2), 75, 'md', 'filled');
        scatter(Covariation.J31(:,1),Covariation.J31(:,2), 75, 'cv', 'filled');
        legend('J12', 'L2', 'J23', 'L3', 'J31')
        title('Map 9mers onto riboswitch')
        box on
        % colorbar
        % xlim([-30, 0])
        % caxis([-1.5, 1.5])
        pbaspect([1 1 1])
        hold off
        saveas(xxx,'covariation_9mers_riboswitch_singlestranded.eps')
    
    end
    
end



%% Figures specific for the ATP RNA Aptamer

xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
hold on
colormap 'cool'
scatter( ...
    reshape(ACEs.deltaG_Tiled5Prime(1:18,1:14),[],1), ...
    reshape(Hyb_data.medianSRS_C_Tiled5Prime(1:18,1:14),[],1), 25, ...
    reshape(ACEs.deltaG_self_Tiled5Prime(1:18,1:14),[],1), 'o', 'filled')
scatter( ...
    reshape(ACEs.deltaG_Tiled5Prime(19:end,1:14),[],1), ...
    reshape(Hyb_data.medianSRS_C_Tiled5Prime(19:end,1:14),[],1), 25, ...
    'ks', 'filled')
caxis([-1.5, 1.5])
title('compare deltaG vs Hyb Signal for 5prime and 3prime ACEs')
box on
legend('5prime', '3prime')
colorbar
xlim([-35, -5])
pbaspect([1 1 1])
hold off
saveas(xxx,'deltaG_vs_KHyb_ATPRNA.eps')




%% Figures Specific for the TBA DNA Aptamer - Sodium Buffer

xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
hold on
colormap 'cool'
scatter( ...
    reshape(ACEs.deltaG_Tiled5Prime(1:21,1:9),[],1), ...
    reshape(Hyb_data.medianSRS_C_Tiled5Prime(1:21,1:9),[],1), 25, ...
    reshape(ACEs.deltaG_self_Tiled5Prime(1:21,1:9),[],1), 'o', 'filled')
scatter( ...
    reshape(ACEs.deltaG_Tiled5Prime(22:end,1:9),[],1), ...
    reshape(Hyb_data.medianSRS_C_Tiled5Prime(22:end,1:9),[],1), 25, ...
    'ks', 'filled')
caxis([-1.5, 1.5])
title('compare deltaG vs Hyb Signal for 5prime and 3prime ACEs for TBA in sodium buffer')
box on
legend('5prime', '3prime')
colorbar
xlim([-25, -5])
pbaspect([1 1 1])
hold off
saveas(xxx,'deltaG_vs_KHyb_TBASodium.eps')





%% Make RED and GREEN channel heat maps for TBA DAs
if CompareRedGreenChannelsFlag == 1;
    
    xxx = figure;
    hold on
    sc(data_REDGREEN.AbsoluteAssay1Tiled5Prime, [0, 5000], 'jet',  [1,1,1])
    title('Absolute signal of thrombin in Sodium')
    colorbar
    hold off
    saveas(xxx,'TBA_REDsignal_absolute_sodium.png')


    xxx = figure;
    hold on
    sc(data_REDGREEN.NormAssay1Tiled5Prime, [0, 1], 'jet',  [1,1,1])
    title('Normalized binding of thrombin in Sodium')
    colorbar
    hold off
    saveas(xxx,'TBA_REDsignal_normalized_sodium.png')

end


%% Compare Sodium and Potassium for TBA DAs.
if CompareSodiumPotassiumFlag == 1;
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    scatter(data_SodiumPotassium.Sodium.Baseline(1:225),data_SodiumPotassium.Potassium.Baseline(1:225), 25, ACEs.deltaG_self(1:225), 'filled')
    title('compare baseline switch rates: Potassium vs. Sodium Bufer')
    box on
    colorbar
    xlim([-0.05, 0.5])
    ylim([-0.05, 0.5])
    caxis([-1.5, 1.5])
    pbaspect([1 1 1])
    hold off
    saveas(xxx,'baseline_Potassium_vs_Sodium_TBA.eps')
    
    
    xxx = figure('DefaultAxesFontName', 'Arial','DefaultAxesFontSize',12);
    hold on
    colormap 'cool'
    scatter(data_SodiumPotassium.Sodium.MaxLigand(1:225),data_SodiumPotassium.Potassium.MaxLigand(1:225), 25, ACEs.deltaG_self(1:225), 'filled')
    title('compare max ligand switch rates: Potassium vs. Sodium Bufer')
    box on
    colorbar
    xlim([-0.05, 0.5])
    ylim([-0.05, 0.5])
    caxis([-1.5, 1.5])
    pbaspect([1 1 1])
    hold off
    saveas(xxx,'maxLigand_Potassium_vs_Sodium_TBA.eps')
    
end


cd ../..

