function simResult = ApplyDFE(simResult,simSettings)
% Applies Decision Feedback Equalizer (DFE) to cancel post-cursor ISI with continuous-time signal correction.
%
% INPUTS:
%   simResult.output                - RX FFE output signal (vector)
%   simResult.eq.thresh             - Decision threshold
%   simResult.eq.tapsDFE            - Vector of post-cursor tap coefficients
%   simSettings.samplesPerSymb      - Number of samples per symbol (sampling rate control)
%
% OUTPUT:
%   outputDFE                       - Continuous-time DFE output signal
    
    % Backup current output
    simResult.outputPrevious = simResult.output;

    output = simResult.output;
    thresh = simResult.eq.thresh;
    tapsDFE = simResult.eq.tapsDFE;
    samplesPerSymb = simSettings.samplesPerSymb;

    % Backup current output
    simResult.outputPrevious = simResult.output;

    % Initialize output signal
    outputDFE = output;

    % Initialize detected symbol array
    detected_symbols = zeros(size(output));
    
    % Align the sampling point with the main cursor
    [~,mainIdx] = max(output);
    n_start = mod(mainIdx, samplesPerSymb);

    % Ensure valid index (at least 1)
    if n_start == 0
        n_start = samplesPerSymb;
    end
    
    % DFE Iterative Calculation
    for n = n_start:samplesPerSymb:length(output)
        % Decision device with threshold
        if outputDFE(n) > thresh
            detected_symbols(n) = 1;
            % disp(n);   % This should be the main cursor
        else
            detected_symbols(n) = -1;
        end

        % DFE in discrete time
        for k = 1:length(tapsDFE)
            delayed_idx = n + (k-1/2) * samplesPerSymb;  % Delay by k*UI
            if delayed_idx <= length(output) && detected_symbols(n) == 1
                % Spread correction effect over symbol interval
                idx_range = delayed_idx : min(delayed_idx + samplesPerSymb - 1, length(output));
                outputDFE(idx_range) = outputDFE(idx_range) - tapsDFE(k) * detected_symbols(n);
            end
        end
    end
    simResult.output = outputDFE;
end
