function [Tiled5Prime, Tiled3Prime, Smoothed, CountsMatrix] = tileSmoothHeatMapData(compiledData, compiledSequences, plotMatrix, startingTilingInAptamer, aptamerSequenceLength, misMatchIndex)
%% TILESMOOTHHEATMATDATA Tile and smooth a dataset in order to make heatmaps
%
%   Required inputs:    1. Dataset to be tiled/smoothed
%                       2. Seqeunces to be tiled/smoothed
%                       3. Matrix outlining how to configure the heatmap
%                          using indexes of mer matches
%                       4. mer length to start tiling
%                       5. Lenght of the aptamer sequence
%                       6. Matrix outlining how to configure the heatmap
%                          using indexes of mer mismatches
%
%   Outputs:            1. Data tiled in a 5' heat map
%                       2. Data tiled in a 3' heat map
%                       3. Data smoothed by linearly projecting into a heat map (if this section of the code is uncommented)
%                       4. Matrix of the nubmer of ACEs making up each tile in the smoothed heat map
%


%% Initialize the vectors
Tiled5Prime(1:aptamerSequenceLength,1:size(plotMatrix,1)) = NaN;
Tiled3Prime(1:aptamerSequenceLength,1:size(plotMatrix,1)) = NaN;
Smoothed(1:aptamerSequenceLength,1:size(plotMatrix,1)) = 0;
CountsMatrix(1:aptamerSequenceLength,1:size(plotMatrix,1)) = 0;

skip = 0;

%% Tile the data into 5' or 3' enantio heat maps
for i = 1:size(plotMatrix,1)
    
    indexes = plotMatrix{i}(startingTilingInAptamer:end);
    
    
    % Control heap map look for mismatches only:
    if plotMatrix{i}(1) == misMatchIndex(1)
        skip = 0;
        
    elseif (plotMatrix{i}(1) > misMatchIndex(1)) && (plotMatrix{i}(1) < misMatchIndex(2))
        skip = skip+1;
        
    elseif plotMatrix{i}(1) == misMatchIndex(3)
        skip = 0;
        
    elseif (plotMatrix{i}(1) > misMatchIndex(3)) && (plotMatrix{i}(1) < misMatchIndex(4))
        skip = skip+1;
        
    elseif plotMatrix{i}(1) >= misMatchIndex(end)
        skip = 0;
    end
    
    
    Tiled5Prime(startingTilingInAptamer + skip : ...
        startingTilingInAptamer + size(indexes,2) - 1 + skip,i) = compiledData(indexes);
    
    Tiled3Prime(end-size(indexes,2)+1 - skip : end - skip,i)= compiledData(indexes);
    
    
    
    
%%     Uncomment for Smoothing (This section of code has been tested on ATP DNA DAs).
    % Smooth the data onto a heatmap of the aptamer sequence vs probe length
% %     for j=1:size(indexes,2)
% %         
% %         first = j + startingTilingInAptamer-1;
% %         
% %         last = first + size(compiledSequences{indexes(j)},2)-1;
% % % 
% % %         if strncmpi(Probes.sequences(i),'TTTTT',5) == 1 || ...
% % %            strncmpi(fliplr(Probes.sequences(i)),'TTTTT',5) == 1
% % %             last = last-10;
% % %         end
% %         
% %         Smoothed5Prime(first:last,i) = Smoothed5Prime(first:last,i) + nansum([0 compiledData(indexes(j))]);
% %         CountsMatrix(first:last,i) = CountsMatrix(first:last,i) + 1;
% %         
% %     end

end

Smoothed(Smoothed==0) = NaN;
Smoothed = Smoothed./CountsMatrix;

end

