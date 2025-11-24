function ffe_taps = AdaptTransmitterFFE(simResult, simSettings)
%   Adapts transmitter FFE taps to cancel precursor interference
%
%   Implements the method in:
%   "High-Speed Wireline Links—Part II: Optimization and Performance Assessment"
%   Section III.C.1 TX FFE Optimization
%
%   INPUTS:
%       simResult.cursorValues      - Array of channel pulse response values
%       simSettings.opt.numTXFFE    - Number of TX FFE taps (N)
%
%   OUTPUT:
%       ffe_taps                    - Optimized and normalized FFE tap coefficients

    cursorValues = simResult.cursorValues;
    numTXFFE = simSettings.opt.numTXFFE;

    % Step 1: Extract segment with main + (N-1) pre-cursors
    [~, mainIdx] = max(cursorValues);
    firstIdx = max(1, mainIdx - numTXFFE + 1);
    h = cursorValues(firstIdx:mainIdx); % h = [h_{-N+1}, ..., h0]
    h = h(:).';  % Ensure row vector

    % Step 2: Construct pulse response matrix H (N × (2N-1))
    N = numTXFFE;
    H = zeros(N, 2*N - 1);
    for i = 1:N
        H(i, i:i + N - 1) = h;
    end

    % Step 3: Construct target pulse response g*
    g_star = zeros(1, size(H, 2));
    g_star(end) = h(end);  % Desired main cursor value at the end

    % Step 4: Compute optimal FFE taps
    % v_opt = g_star * H^T * inv(H * H^T)
    ffe_taps = g_star * H' / (H * H');

    % Step 5: Normalize to maintain consistent TX swing
    ffe_taps = ffe_taps / sum(abs(ffe_taps));
end

% function ffe_taps = AdaptTransmitterFFE(simResult, simSettings)
% % GPU-accelerated TX FFE tap adaptation
% 
%     cursorValues = gpuArray(simResult.cursorValues);  % ensure GPU
%     numTXFFE = simSettings.opt.numTXFFE;
% 
%     % Step 1: Extract [h_{-N+1}, ..., h_0]
%     [~, mainIdx] = max(cursorValues);
%     firstIdx = max(1, mainIdx - numTXFFE + 1);
%     h = cursorValues(firstIdx:mainIdx);  % 1D gpuArray
%     h = h(:).';  % Ensure row vector
% 
%     % Step 2: Construct H (N × (2N−1)) from h
%     N = numTXFFE;
%     H = gpuArray.zeros(N, 2 * N - 1);
%     for i = 1:N
%         H(i, i:i + N - 1) = h;
%     end
% 
%     % Step 3: Construct target pulse response g_star
%     g_star = gpuArray.zeros(1, size(H, 2));
%     g_star(end) = h(end);  % target = [0 ... 0 h_0]
% 
%     % Step 4: Solve for optimal taps
%     ffe_taps = g_star * H' / (H * H');
% 
%     % Step 5: Normalize
%     ffe_taps = ffe_taps / sum(abs(ffe_taps));
% end