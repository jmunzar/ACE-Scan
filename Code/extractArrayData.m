function [Hyb_rawdata, Assay_rawdata] = extractArrayData(filename_Hyb, filename_Assay, workingfolder, N_SubArrays)
%% EXTRACTARRAYDATA Extract the data in the ACE scanning microarray .gal file
%
%   Required inputs:    1. Hyb intensity file (.csv)
%                       2. Assay intensity file (.csv)
%                       3. Working folder
%                       4. Number of subarrays in the microarray
%
%   Outputs:            1. Hyb array raw data, structured
%                       2. Assay array raw data, structured


%%  Delete later
% % clc
% % clear all
% % close all
% % filename_Hyb = '2016.02.25_Slide2_10.10_ATP_AfterHyb_G_gain5_2micron_SLOT01_S01 transposed rotated 180_30pixel.csv';
% % filename_Assay = '2016.02.25_Slide2_10.10_ATP_AfterAssay_G_gain20_XDR_20_2micron_SLOT01_S01_L transposed rotated 180_30pixel.csv';
% % workingfolder = '/Users/jmunzar/DJGROUP/Aptamers/Software/ATP_10.10_Results';
% % N_SubArrays = 6;


%%  Read in files
StartRow = 1;
StartCol = 2;
% Open files, read file contents
try
    Hyb_Array_AllData = csvread(filename_Hyb,StartRow,StartCol);
catch me
    disp('Can not find Hyb file')
end
try
    Assay_Array_AllData = csvread(filename_Assay,StartRow,StartCol);
catch me
    disp('Can not find Assay file')
end


%%  Sort long Array dataset into an ordered array (i = ACE number, j = subarray)
iIndex = 0;
jIndex = 1;
for i = 1:size(Hyb_Array_AllData,1)
    iIndex = iIndex + 1; % count iIndex as we go
    if iIndex == size(Hyb_Array_AllData,1)/N_SubArrays+1
        iIndex = 1; %reset iIndex counter when a full SubArray is done
        jIndex = jIndex+1;% and increment the jIndex
    end
    
    Temp1h(iIndex,jIndex) = Hyb_Array_AllData(i,1);
    Temp2h(iIndex,jIndex) = Hyb_Array_AllData(i,2);
    Temp4h(iIndex,jIndex) = Hyb_Array_AllData(i,4);
    Temp5h(iIndex,jIndex) = Hyb_Array_AllData(i,5);
    Temp7h(iIndex,jIndex) = Hyb_Array_AllData(i,7);
    Temp8h(iIndex,jIndex) = Hyb_Array_AllData(i,8);
    
    Temp1a(iIndex,jIndex) = Assay_Array_AllData(i,1);
    Temp2a(iIndex,jIndex) = Assay_Array_AllData(i,2);
    Temp4a(iIndex,jIndex) = Assay_Array_AllData(i,4);
    Temp5a(iIndex,jIndex) = Assay_Array_AllData(i,5);
    Temp7a(iIndex,jIndex) = Assay_Array_AllData(i,7);
    Temp8a(iIndex,jIndex) = Assay_Array_AllData(i,8);
    
end

%%  Store all data into two structures for Hyb and Assay microarrays

Hyb_rawdata = struct(   'meanSignal', Temp1h, ...
                        'meanBackground', Temp2h, ...
                        'medianSignal', Temp4h, ...
                        'medianBackground', Temp5h, ...
                        'stdevSignal', Temp7h, ...
                        'stdevBackground', Temp8h);

Assay_rawdata = struct( 'meanSignal', Temp1a, ...
                        'meanBackground', Temp2a, ...
                        'medianSignal', Temp4a, ...
                        'medianBackground', Temp5a, ...
                        'stdevSignal', Temp7a, ...
                        'stdevBackground', Temp8a);

end