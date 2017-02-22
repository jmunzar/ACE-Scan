function [TiledSequences, Tm_TiledSequences5, DeltaG, DeltaG_self] = extractOligoProperties(filename_ACE_Sequences, N_ACEs) 
%% EXTRACTOLIGOPROPERTIES Extract Tm data from the .csv list of ACE Sequences
%
%  Required inputs:     1. .csv with the complete list of ACE sequences
%                           in 5' to 3' format, also contains deltaG data
%                       2. Number of ACEs in the file
%                       
%   Outputs:            1. ACE sequences in an array
%                       2. Tm of each ACE accoring to Nearest Neighbor (NN) matlab model
%                       3. deltaG value for each ACE hybridizing to the
%                          aptamer accoring to UNAfold (TwoState application).
%                       4. deltaG value for each ACE hybridizing to
%                          itself using UNAfold (Zipfold application).



% Open the file, copy file contents
try
    fid=fopen(filename_ACE_Sequences,'r');
    if fid == -1
        error ('LS:NoImageFile','no gal .csv file present')
    end
str2=textscan(fid,'%s %f %f', 'Delimiter', ',');
catch me
    disp('Can not read file')
end
if fid ~= -1;
    fclose(fid);
end

% Find Tm and deltaG from .csv or predicted by matlab.
TiledSequences=cellstr(strvcat(str2{1}));
DeltaG = str2{2};
DeltaG_self = str2{3};

for i=1:N_ACEs
    N_basecount = basecount(TiledSequences{i}); % Dont return
    GCcontent(i,1) = (N_basecount.G + N_basecount.C) / (N_basecount.A + N_basecount.T + N_basecount.G + N_basecount.C); % Dont return
    Length_TiledSequences(i,1) = N_basecount.A + N_basecount.T + N_basecount.G + N_basecount.C; % Dont return
    Tm = oligoprop(TiledSequences{i},'Salt',0.6,'Primerconc',1e-6);
%     Tm_TiledSequences1(i,1) = Tm.Tm(1); Other models of DNA hybridization
%     Tm_TiledSequences2(i,1) = Tm.Tm(2);
%     Tm_TiledSequences3(i,1) = Tm.Tm(3);
    Tm_TiledSequences5(i,1) = Tm.Tm(5);
end


end