function [ Km, Km_Std, Vmax, Vmax_Std, MSError] = extractMichaelisMenten( EnzymeArray, LigandConc, Vmax_min, Vmax_max )
%% EXTRACT MICHAELISMENTEN Extract Km and Vmax for an enzyme reaction
%
%  Required inputs:     1. 2D array of enzyme velocities, each column corresponding to a concentration
%                       2. Vector of corresponding ligand concentratons
%                       3. Cutoff for minimum Vmax value to return
%                       4. Cutoff for maximum Vmax value to return
%
%   Outputs:            1. Km, the michaelis constant
%                       2. Km_Std, 95% CI for Km
%                       3. Vmax, the maximum turnover rate at infinite ligand concentraton
%                       4. Vmax_Std, 95% CI for Vmax
%                       5. MSError, mean square error of the fitted model
%
%   This function is adapted from 
%   https://sakai.duke.edu/access/content/group/25e08a3d-9fc4-41b0-a7e9-815732c1c4ba/New%20folder/Stat%20Topic%20Files/Non-Linear%20Regression/EnzymeKineticsExample.pdf


%% Assign initial data vectors
Km = zeros(1,(size(EnzymeArray,1)));
Vmax = zeros(1,(size(EnzymeArray,1)));

%% Extract values for each ACE
for ACE = 1:size(EnzymeArray,1)
    
    % Define Vo(koff) using buffer only rates.
    vo=-(EnzymeArray(ACE,1:end-1)-EnzymeArray(ACE,end)); % Effective enzyme velocity

    
    % Calculate unknown coefficients in the model using nlinfit with 2 free
    % parameters (Vmax, KFit).
    model=@(b,x) b(1).*LigandConc./(b(2)+LigandConc);
    
    % Apply initial guess. If many errors are reported, consider trying
    % different initial guesses.
    initialguess=[-(EnzymeArray(ACE,1) - EnzymeArray(ACE,end)) .0001];
    
    
    try
        [betaB,R,J,CovB,MSE] = nlinfit(LigandConc,vo,model,initialguess);
        betaB;
        % 95% CI of coeffieicents
        betaBci = nlparci(betaB,R,J);
        
    catch
        fail(ACE) = 1;
    end
    
    % Put constraints on the values of Vmax, Km that can be accepted
    % (ONly report good fits with Vmax between defined min and max and Km between min and max ligand
    % concentrations)
    
    % Also calculate 95% confidence intervals on variables
        
    if exist('betaB') == 0
        exist('betaB')
        MSError(ACE) = NaN;
        Vmax(ACE) = NaN;
        Vmax_Std(ACE) = NaN;
        Km(ACE) = NaN;
        Km_Std(ACE) = NaN;
        
    elseif betaB(1) > Vmax_min && betaB(1) < Vmax_max ...
            && betaB(2)>LigandConc(end) &&  betaB(2)<LigandConc(1)
        MSError(ACE) = MSE;
        Vmax(ACE) = betaB(1);
        Vmax_Std(ACE) = (betaBci(3) - betaB(1))/2;
        Km(ACE) = betaB(2);
        Km_Std(ACE) = (betaBci(4) - betaB(2))/2;
        
    else
        MSError(ACE) = NaN;
        Vmax(ACE) = NaN;
        Vmax_Std(ACE) = NaN;
        Km(ACE) = NaN;
        Km_Std(ACE) = NaN;
    end
    
end

MSError = MSError';
Vmax = Vmax';
Vmax_Std = Vmax_Std';
Km = Km';
Km_Std = Km_Std';

end