function simResult = ApplyFFE(simResult, tapsFFE, simSettings)
% Applies the optimized FFE taps to the main signal.
%
% INPUTS:
%   simResult.output              - Main channel response signal
%   tapsFFE                       - Optimized FFE tap coefficients
%   simSettings.samplesPerSymb    - Number of samples per symbol (UI)
%
% OUTPUT:
%   simResult                     - Updated simResult with FFE-applied signals

    % Backup current output
    simResult.outputPrevious = simResult.output;

    % Upsample FFE taps to match oversampling rate
    samplesPerSymb = simSettings.samplesPerSymb;
    tapsFFE_upsampled = zeros(1, (length(tapsFFE) - 1) * samplesPerSymb + 1);
    tapsFFE_upsampled(1:samplesPerSymb:end) = tapsFFE;

    % Apply FFE to main signal
    simResult.output = conv(simResult.output, tapsFFE_upsampled, 'same');
end

% function simResult = ApplyFFE(simResult, tapsFFE, simSettings)
% % GPU-accelerated FFE application to oversampled waveform
% 
%     % Backup current output
%     simResult.outputPrevious = simResult.output;
% 
%     % Ensure GPU input
%     output = simResult.output;                     % already on GPU
%     tapsFFE = gpuArray(tapsFFE);                   % convert taps to GPU
%     samplesPerSymb = simSettings.samplesPerSymb;
% 
%     % Upsample taps to match oversampling rate
%     upLen = (length(tapsFFE) - 1) * samplesPerSymb + 1;
%     tapsFFE_upsampled = gpuArray.zeros(1, upLen);
%     tapsFFE_upsampled(1:samplesPerSymb:end) = tapsFFE;
% 
%     % Apply convolution on GPU
%     simResult.output = conv(output, tapsFFE_upsampled, 'same');
% end
