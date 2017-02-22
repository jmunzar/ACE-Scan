%MYcroArray_Generator.m
%
%

% Clear matlab variables.
clc;    
clear all;
close all;
workspace;
pwd;



%% 1. SET THE DIRECTORY - update the following line to point to the location of this folder on your system
cd('/Users/jmunzar/GitHub/ACE-Scanning/Code/');




%% 2. Define Aptamer sequence

% RNA ATP aptamer as an example, with added 3 bp at 5' and 3' ends taken from primers used in the original SELEX experiment:
% CTT --- GGGTTGGGAAGAAACTGTGGCACTTCGGTGCCAGCAACCC --- CAT
% Input as a DNA sequence:
AptamerVanillaSequence = 'CTTGGGTTGGGAAGAAACTGTGGCACTTCGGTGCCAGCAACCCCAT'   %Aptamer sequence, in 5' -> 3' order




%% 3. Define ACE designs to be tested
%%% Note that individual subsections of the code below can be commented if
%%% the design of ACEs isn't required.

StartingMerNoMisMatch = 7; % Starting size of N-mer ACEs
EndingMerNoMisMatch = 20; % Ending size of N-mer ACEs (for all, set to size(AptamerVanillaSequence,2);
MerSingleMisMatchArray = [12 15]; %% sizes of N-mer ACEs with 1 mismatch






%% Seed random nubmer generator
s = RandStream('mt19937ar','Seed',1); %Use the same SEED for the random nubmer generator for each DA family.
RandStream.setGlobalStream(s);




%% 4. Generate all N-mers of perfect overlap

k = 0;
for Mer = StartingMerNoMisMatch:EndingMerNoMisMatch %Loop from the smallest to the largest length of the ACE
    LengthLoop_EachNmer = size(AptamerVanillaSequence,2) - Mer + 1;
    
    
    for i = 1:LengthLoop_EachNmer %Loop for each N-mer size of ACE
        
        
        ArraySequences_Master(k+i,1:Mer) = AptamerVanillaSequence(i:i+Mer-1);     % Take partial sequence data from the aptamer (tiling 5'->3')
        
        
    end
    % ArraySequences_Master
    k = k+LengthLoop_EachNmer;
    
end



%% 5. Generate all N-mers of perfect overlap with poly A extensions % % %
%%% Comment this section to remove these ACEs in the design

%%% Adds a Poly 10 A to 3' or to 5' ends, which will be turned into a polyT at
%%% 5' or 3' ends after inversing the ACE sequences.



AddOnPolyA = 'AAAAAAAAAA';
AddOnSize = size(AddOnPolyA,2);

%%% Add 5' polyT
for Mer = StartingMerNoMisMatch:EndingMerNoMisMatch %Loop from the smallest to the largest length of the ACE
    LengthLoop_EachNmer = size(AptamerVanillaSequence,2) - Mer + 1;
    
    for i = 1:LengthLoop_EachNmer %Loop for each N-mer size of ACE
        
        ArraySequences_Master(k+i,1:Mer+AddOnSize) = strcat(AptamerVanillaSequence(i:i+Mer-1),AddOnPolyA); % Take partial sequence data from the aptamer (tiling 5'->3')
        
    end
    % ArraySequences_Master
    k = k+LengthLoop_EachNmer;
    
end

%%% Add 3' polyY
for Mer = StartingMerNoMisMatch:EndingMerNoMisMatch %Loop from the smallest to the largest length of the ACE
    LengthLoop_EachNmer = size(AptamerVanillaSequence,2) - Mer + 1;
    
    for i = 1:LengthLoop_EachNmer %Loop for each N-mer size of ACE
        
        ArraySequences_Master(k+i,1:Mer+AddOnSize) = strcat(AddOnPolyA, AptamerVanillaSequence(i:i+Mer-1)); % Take partial sequence data from the aptamer (tiling 5'->3')
        
    end
    % ArraySequences_Master
    k = k+LengthLoop_EachNmer;
    
end




%% 6. Generate all N-mers with 1-base mismatch

for Mer = MerSingleMisMatchArray %Loop from the smallest to the largest length of the ACE
    
    LengthLoop_EachNmer = (size(AptamerVanillaSequence,2) - Mer + 1);
    
    for i = 1:LengthLoop_EachNmer %Loop for each N-mer size of ACE
        
        for j = 1:Mer
            k = k+1;
            ArraySequences_Master(k,1:Mer) = AptamerVanillaSequence(i:i+Mer-1);  % Take partial sequence data from the aptamer (tiling 5'->3')
            
            if ArraySequences_Master(k,j) == 'T'
                ArraySequences_Master(k,j) = 'A';
            else
                ArraySequences_Master(k,j) = 'T';
            end

        end
        
    end

end





%% 7. Add long consensus ACEs as P.C. sequences

for Mer = size(AptamerVanillaSequence,2)-5:size(AptamerVanillaSequence,2) %Loop from the smallest to the largest length of the ACE
    LengthLoop_EachNmer = size(AptamerVanillaSequence,2) - Mer + 1;
    
    for i = 1:LengthLoop_EachNmer %Loop for each N-mer size of ACE
        
        ArraySequences_Master(k+i,1:Mer) = AptamerVanillaSequence(i:i+Mer-1);     % Take partial sequence data from the aptamer (tiling 5'->3')
                
    end

    k = k+LengthLoop_EachNmer;
    
end





%% 8. Take reverse complement of each blocker


ArraySequences_Final = ArraySequences_Master;
for i = 1:size(ArraySequences_Master,1)
    

    ArraySequences_Final(i,:) = regexprep(ArraySequences_Final(i,:),'A','B'); % Replace A, C, G, T's with T, G, C, A's
    ArraySequences_Final(i,:) = regexprep(ArraySequences_Final(i,:),'T','A');
    ArraySequences_Final(i,:) = regexprep(ArraySequences_Final(i,:),'B','T');
    ArraySequences_Final(i,:) = regexprep(ArraySequences_Final(i,:),'C','D');
    ArraySequences_Final(i,:) = regexprep(ArraySequences_Final(i,:),'G','C');
    ArraySequences_Final(i,:) = regexprep(ArraySequences_Final(i,:),'D','G');
    
    ArraySequences_Final(i,:) = fliplr(ArraySequences_Final(i,:));            % Reverse the sequence into 5' -> 3' order
    
end



%% 9. Shuffle and UnShuffle the array of ACEs


Shuffle = randperm(size(ArraySequences_Final,1));

ArraySequences_Shuffled = ArraySequences_Final(Shuffle,:);


UnShuffle = zeros(size(Shuffle,1)); %
for i = 1:size(Shuffle,2)
    UnShuffle(Shuffle(i)) =  i;
end

ArraySequences_UnShuffled = ArraySequences_Shuffled(UnShuffle,:);




%% 10. Write all ACEs to files

S = cellstr(ArraySequences_Shuffled);
S2 = cellstr(ArraySequences_UnShuffled);


SS = strtrim(S);  %remove whitespace in the array
SS2 = strtrim(S2);  %remove whitespace in the array


%%% Uncomment this break to stop before writing over the files
% break

mkdir('ACEdesigns')
cd('ACEdesigns')

dlmwrite('OligoListShuffled.csv',char(SS),'delimiter','');
dlmwrite('OligoListUnshuffled.csv',char(SS2),'delimiter','');


xlswrite('Test.xls',char(SS))
xlswrite('Test2.xls',char(SS2))

fid = fopen('Shuffled_Final.txt','wt');
   for i = 1 : size(SS,1)
        fprintf(fid,'%s\n',SS{i,:});
   end
fclose(fid);


fid = fopen('UnShuffled_Final.txt','wt');
   for i = 1 : size(SS2,1)
        fprintf(fid,'%s\n',SS2{i,:});
   end
fclose(fid);


save('Database');
cd ..
