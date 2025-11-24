function [tapsFFE, optSamplingPoint] = AdaptReceiverFFE(simResult, simSettings)
% Computes optimized RX FFE taps and optimal sampling point.
%
% INPUTS:
%   simResult.cursorValues          - Response cursor values
%   simResult.eq.tapsDFE            - Post-tap factor for DFE co-optimization
%   simSettings.opt.numRXFFE        - Number of RX FFE taps
%
% OUTPUTS:
%   tapsFFE                         - Optimized RX FFE tap coefficients
%   optSamplingPoint                - Optimal sampling point index

    % Compute RX noise power shaped by RX BW limitation, CTLE, and VGA
    rxNoisePower = ComputeRxNoisePower(simResult, simSettings);
    % disp(rxNoisePower);

    % Extract settings and inputs
    cursorValues = simResult.cursorValues;
    tapsDFE = simResult.eq.tapsDFE;
    numRXFFE = simSettings.opt.numRXFFE;

    % Initialization
    minMMSE = Inf;
    optSamplingPoint = 0;
    w_opt = zeros(numRXFFE, 1);

    % Construct matrix X (Toeplitz-like)
    X = zeros(numRXFFE, numRXFFE + length(cursorValues) - 1);
    for row = 1:numRXFFE
        X(row, row : row + length(cursorValues) - 1) = cursorValues;
    end
    y_length = size(X, 2);

    % Sweep over candidate sampling points
    for i = 1:y_length
        % Build ideal target response
        y = zeros(y_length, 1);
        y(i) = 1;

        % Add post-cursors for DFE co-optimization
        for j = 1:length(tapsDFE)
            idx = i + j;
            if idx <= y_length
                y(idx) = tapsDFE(j);
            end
        end

        % Build system A w = b
        A = X * X.';
        b = X * y;

        % Check conditioning and apply regularization if needed
        if rcond(A) < 1e-6
            reg_strength = rxNoisePower;
            A = A + reg_strength * eye(size(A));
        end

        % Solve system
        try
            w_i = A \ b;
        catch
            continue;  % Skip this sampling point if solving fails
        end

        % Compute MMSE
        y_hat = X.' * w_i;
        mmse = norm(y_hat - y)^2 + rxNoisePower * norm(w_i)^2;

        % Track best result
        if mmse < minMMSE
            minMMSE = mmse;
            w_opt = w_i;
            optSamplingPoint = i;
        end
    end

    % Normalize
    if norm(w_opt) > 0
        tapsFFE = w_opt / norm(w_opt);
    else
        tapsFFE = w_opt;
    end
end

% function [tapsFFE, optSamplingPoint] = AdaptReceiverFFE(simResult, simSettings)
% % GPU-accelerated RX FFE adaptation with MMSE + DFE co-optimization
% 
%     % Compute RX noise power (remains CPU-side, since FFT-based)
%     rxNoisePower = ComputeRxNoisePower(simResult, simSettings);
% 
%     cursorValues = gpuArray(simResult.cursorValues);
%     tapsDFE = gpuArray(simResult.eq.tapsDFE);
%     numRXFFE = simSettings.opt.numRXFFE;
% 
%     minMMSE = Inf;
%     optSamplingPoint = 0;
%     w_opt = gpuArray.zeros(numRXFFE, 1);
% 
%     % Construct matrix X on GPU (Toeplitz-like rows of cursorValues)
%     X = gpuArray.zeros(numRXFFE, numRXFFE + length(cursorValues) - 1);
%     for row = 1:numRXFFE
%         X(row, row : row + length(cursorValues) - 1) = cursorValues;
%     end
%     y_length = size(X, 2);
% 
%     for i = 1:y_length
%         % Target pulse response
%         y = gpuArray.zeros(y_length, 1);
%         y(i) = 1;
% 
%         % Add DFE postcursors
%         for j = 1:length(tapsDFE)
%             idx = i + j;
%             if idx <= y_length
%                 y(idx) = tapsDFE(j);
%             end
%         end
% 
%         % Solve A w = b
%         A = X * X';
%         b = X * y;
% 
%         if rcond(gather(A)) < 1e-6
%             A = A + rxNoisePower * eye(size(A), 'gpuArray');
%         end
% 
%         try
%             w_i = A \ b;
%         catch
%             continue;
%         end
% 
%         y_hat = X' * w_i;
%         mmse = norm(y_hat - y)^2 + rxNoisePower * norm(w_i)^2;
% 
%         if mmse < minMMSE
%             minMMSE = mmse;
%             w_opt = w_i;
%             optSamplingPoint = i;
%         end
%     end
% 
%     % Normalize & return to CPU
%     if norm(w_opt) > 0
%         tapsFFE = gather(w_opt / norm(w_opt));
%     else
%         tapsFFE = gather(w_opt);
%     end
% end

function sigma_n2 = ComputeRxNoisePower(simResult, simSettings)
% Computes the noise power σ_n² at the input of the RX FFE by shaping RX white noise
% through receiver bandwidth, CTLE, and VGA frequency responses.
%
% INPUTS:
%   simSettings.rxPSD            - RX noise power spectral density (V^2/Hz)
%   simResult.H_bw               - RX BW Limitation Transfer function
%   simResult.H_ctle             - CTLE Transfer function
%   simResult.H_vga              - VGA Transfer function
%
% OUTPUT:
%   sigma_n2                     - Noise power at RX FFE input

    % === Basic Parameters ===
    rxPSD = simSettings.rxPSD * 1e-9;
    df = simResult.fftFreqs(2) - simResult.fftFreqs(1);

    % === 1. RX BW Transfer Function ===
    H_bw = simResult.H_bw;

    % === 2. CTLE Transfer Function ===
    H_ctle = simResult.H_ctle;

    % === 3. VGA Transfer Function ===
    H_vga = simResult.H_vga;

    % === Total Frequency Response ===
    H_total = abs(H_bw .* H_ctle .* H_vga).^2;

    % === Noise Power Integration ===
    sigma_n2 = rxPSD * sum(H_total) * df;
end

