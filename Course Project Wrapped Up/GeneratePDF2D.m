function simResult = GeneratePDF2D(simResult, simSettings, pltSettings, simResources)
% Generates and displays a 2-D PDF over amplitude
% Inputs:
%   - simResult         : Simulation result includes output
%   - simSettings       : Simuation settings
%   - pltSettings       : Plot settings updated from simSettings
%   - simResources      : Struct with crosstalk rational models (rS11, rS12, rS21, rS22)
%
% Output:
%   - simResult         : Updated simulation result with all the 2D pdfs

    % Generate and plot 2-D PDFs
    % Generate ISI pdf
    simResult = GeneratePdfISI(simResult, simSettings, pltSettings);
    
    % Generate Crosstalk pdf
    simResult = GeneratePdfCrosstalk(simResult, simSettings, pltSettings, simResources);
    
    % Generate TX and RX Noise pdf
    simResult = GeneratePdfNoise(simResult, simSettings, pltSettings);
    
    % Generate Jitter pdf
    simResult = GeneratePdfJitter(simResult, simSettings, pltSettings);
    
    % Generate Total pdf
    simResult = GeneratePdfTotal(simResult, simSettings, pltSettings, true);
end

