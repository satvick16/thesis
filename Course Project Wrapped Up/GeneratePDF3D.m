function pdf3D = GeneratePDF3D(simResult, simSettings, pltSettings, simResources)
% Generates a 3D PDF over time and amplitude.
% Inputs:
%   - simResult         : Simulation result includes output
%   - simSettings       : Simuation settings
%   - pltSettings       : Plot settings updated from simSettings
%   - simResources      : Struct with crosstalk rational models (rS11, rS12, rS21, rS22)

% Output:
%   - pdf3D             : Generated 3D pdf plot

    output = simResult.output;
    samplesPerSymb = simSettings.samplesPerSymb;    
    trackPre = simSettings.trackPre;
    trackPost = simSettings.trackPost;
    pdfRange = pltSettings.pdfRange;
    dV = pltSettings.pdfdx;
    dT = 1 / samplesPerSymb;

    % Initialize the 3D PDF matrix
    pdf3D = zeros(samplesPerSymb, length(pdfRange));

    % Find the index of the main cursor (maximum value in dfe_output)
    [~, mainIdx] = max(output);

    % Generate the indices for all cursors
    cursorIndices = (mainIdx - trackPre*samplesPerSymb : samplesPerSymb : mainIdx + trackPost*samplesPerSymb);

    % Ensure we stay within bounds of dfe_output
    cursorIndices = cursorIndices(cursorIndices >= 1 & cursorIndices <= length(output));

    % Time shifts corresponding to oversampling points within one UI
    shift_indices = -floor(samplesPerSymb/2):floor(samplesPerSymb/2);
    
    % Generate Crosstalk pdf
    simResult = GeneratePdfCrosstalk(simResult,simSettings,pltSettings,simResources);
    
    % Generate TX and RX Noise pdf
    simResult = GeneratePdfNoise(simResult,simSettings,pltSettings);

    % === parfor over oversampling phases ===
    parfor t_idx = 1:samplesPerSymb
        % Create local copies of mutable structures
        localSimResult = simResult;
    
        % Shift
        shift = shift_indices(t_idx);
        shiftedIndices = cursorIndices + shift;
        shiftedIndices = shiftedIndices(shiftedIndices >= 1 & shiftedIndices <= length(output));
        localSimResult.cursorValues = output(shiftedIndices);
    
        % ISI, jitter, total pdf
        localSimResult = GeneratePdfISI(localSimResult, simSettings, pltSettings);
        localSimResult = GeneratePdfJitter(localSimResult, simSettings, pltSettings);
        localSimResult = GeneratePdfTotal(localSimResult, simSettings, pltSettings, false);
    
        % Store pdf slice
        localPdfSlice = localSimResult.pdfTotal;
    
        % Store outputs
        pdf3D(t_idx, :) = localPdfSlice;
    end

    total_mass = sum(pdf3D, 'all') * dV * dT;   % Compute full integral
    pdf3D = pdf3D / total_mass;                 % Normalize entire surface
end