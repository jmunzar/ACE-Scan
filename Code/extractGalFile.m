function [N_Rows, N_Columns, N_SubArrays, N_ACEs, N_ReplicatesPerACE, SpotIDnumber_Left, SpotIDnumber_Right, sort_indices_Left, sort_indices_Right] = extractGalFile(filename_Gal, workingfolder) 
%% EXTRACTGALFILE Extract the spot informtation stored in the microarray .gal file
%
%   Required inputs:    1. MYcroarray file (.csv equivalent of the .gal from MYcroarray)
%                       2. Working folder
%
%   Outputs:            1. Number of rows in the microarray
%                       2. Number of columns in the microarray
%                       3. Number of subarrays in the microarray
%                       4. Number of unique ACEs synthesized
%                       5. Number of replicates of each unique ACE
%                       6. The list of Spot ID's on LHS of the microarray
%                       7. The list of Spot ID's on RHS of the microarray
%                       8. A matrix to sort the ACEs accoring to ID for
%                          LHS of the microarray
%                       9. A matrix to sort the ACEs accoring to ID for
%                          RHS of the microarray



%% Open Files

cd(workingfolder); % Move to new working directory
% Open the file, copy file contents

try
    fid=fopen(filename_Gal,'r');
    if fid == -1
        error ('LS:NoImageFile','no gal .csv file present')
    end
    str_gal=textscan(fid,'%s');
catch me
    disp('Can not read file. Check file name or directory')
end
if fid ~= -1;
    fclose(fid);
end


%% Split up data from gal file into a cell data type
splitData_galFile = regexp(cellstr(strvcat(str_gal{1})),'[\,]*','split');



%% Search for start of position data in .gal file
sizecounter = size(splitData_galFile{end},2); % Check for size of last array element

for i = 1:size(splitData_galFile,1)
   if size(splitData_galFile{i},2) == sizecounter
       startingrowSplitData = i+1;
       break
   end
   if i == size(splitData_galFile,1)
        error ('Error: Could not Locate start of position data in gal file, consider setting startingrowSplitData to 7')
   end
end




%% Look at how many total blocks there are:

N_Rows = str2num(splitData_galFile{end}{2});
N_Columns = str2num(splitData_galFile{end}{3});
N_SubArrays = str2num(splitData_galFile{end}{1});

% Find numSpots:
numSpots = size(splitData_galFile,1)-startingrowSplitData+1;
numSpotsPerArray = numSpots/N_SubArrays;

%   Should have 7448 spots per array, N replicates of NN that are real 
%   sequences, rest are control sequences from MYcroarray



%% Extract .gal spot informtation into lists
% SpotIDnumber_Left = zeros(numSpotsPerArray,1);
% SpotIDnumber_Right = zeros(numSpotsPerArray,1);
for i = 1:numSpotsPerArray
    
    counter_Left  = i+startingrowSplitData-1;
	counter_Right = i+startingrowSplitData-1+numSpotsPerArray;
    
    SpotIDnumber_Left{i,1} = splitData_galFile{counter_Left}{4};
    SpotIDnumber_Right{i,1} = splitData_galFile{counter_Right}{4};
    
end



%% Determine how many replicates of each condition are in each block:
% use !isnan and str2double to clear non-number filled cells
[N_ReplicatesPerSpotIndex, temp]=hist(str2double(SpotIDnumber_Left(~isnan(str2double(SpotIDnumber_Left)))),unique(str2double(SpotIDnumber_Left(~isnan(str2double(SpotIDnumber_Left))))));
N_ACEs = size(N_ReplicatesPerSpotIndex,2);
N_ReplicatesPerACE = mode(N_ReplicatesPerSpotIndex);



%% Rename the non-array elements to enable conversion to double for sorting
for i = 1:size(SpotIDnumber_Left,1)
    if size(SpotIDnumber_Left{i},2) == 12
        if SpotIDnumber_Left{i} == 'Ctrl-Pos-555'
            SpotIDnumber_Left{i} = num2str(N_ACEs+1);
        elseif SpotIDnumber_Left{i} == 'Ctrl-Pos-647'
            SpotIDnumber_Left{i} = num2str(N_ACEs+2);
        end
    end
    if size(SpotIDnumber_Left{i},2) == 18
        if SpotIDnumber_Left{i} == 'Ctrl-Stringent-555'
            SpotIDnumber_Left{i} = num2str(N_ACEs+3);
        elseif SpotIDnumber_Left{i} == 'Ctrl-Stringent-647'
            SpotIDnumber_Left{i} = num2str(N_ACEs+4);
        end
    end
    if size(SpotIDnumber_Left{i},2) == 5
        if SpotIDnumber_Left{i} == 'MY-QC'
             SpotIDnumber_Left{i} = num2str(N_ACEs+5);
        end
    end 
    if size(SpotIDnumber_Left{i},2) == 5
        if SpotIDnumber_Left{i} == 'Empty'
             SpotIDnumber_Left{i} = num2str(N_ACEs+6);
        end
    end
    if size(SpotIDnumber_Left{i},2) == 11
        if SpotIDnumber_Left{i} == 'Empty_NCTRL'
             SpotIDnumber_Left{i} = num2str(N_ACEs+6);
        end
    end
end


%% Sort the right hand of the microarray seperately (Due to different layout in gal file)
for i = 1:size(SpotIDnumber_Right,1)
    if size(SpotIDnumber_Right{i},2) == 12
        if SpotIDnumber_Right{i} == 'Ctrl-Pos-555'
            SpotIDnumber_Right{i} = num2str(N_ACEs+1);
        elseif SpotIDnumber_Right{i} == 'Ctrl-Pos-647'
            SpotIDnumber_Right{i} = num2str(N_ACEs+2);
        end
    end
    if size(SpotIDnumber_Right{i},2) == 18
        if SpotIDnumber_Right{i} == 'Ctrl-Stringent-555'
            SpotIDnumber_Right{i} = num2str(N_ACEs+3);
        elseif SpotIDnumber_Right{i} == 'Ctrl-Stringent-647'
            SpotIDnumber_Right{i} = num2str(N_ACEs+4);
        end
    end
    if size(SpotIDnumber_Right{i},2) == 5
        if SpotIDnumber_Right{i} == 'MY-QC'
             SpotIDnumber_Right{i} = num2str(N_ACEs+5);
        end
    end 
    if size(SpotIDnumber_Right{i},2) == 5
        if SpotIDnumber_Right{i} == 'Empty'
             SpotIDnumber_Right{i} = num2str(N_ACEs+6);
        end
    end
    if size(SpotIDnumber_Right{i},2) == 11
        if SpotIDnumber_Right{i} == 'Empty_NCTRL'
             SpotIDnumber_Right{i} = num2str(N_ACEs+6);
        end
    end
end


%% Convert spot IDs to doubles, sort accordingly
spotID_Doubles_Left = str2double(SpotIDnumber_Left);
[spotID_sorted_Left, sort_indices_Left] = sort(spotID_Doubles_Left);

spotID_Doubles_Right = str2double(SpotIDnumber_Right);
[spotID_sorted_Right, sort_indices_Right] = sort(spotID_Doubles_Right);

 
end