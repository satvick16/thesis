function simResult = ApplyReceiverBW(simResult)
% Applies a 4th-order Butterworth low-pass filter to model the receiver bandwidth limitation.
%
% INPUTS:
%   simResult.Yf                 - Input signal in frequency domain
%   simResult.rxBW               - Cutoff frequency (Hz)
%
% OUTPUT:
%   simResult.H_bw               - Transfer function
%   simResult.Yf                 - Equalized signal in frequency domain
%   simResult.output             - Equalized signal in time domain

    % Backup previous output
    simResult.outputPrevious = simResult.output;

    % Extract parameters
    Yf = simResult.Yf;
    x = simResult.fftFreqs(:) / simResult.rxBW;

    % 4th-order Butterworth
    H_bw = 1 ./ (1 - 3.4142 * x.^2 + x.^4 + 1j * 2.6131 * (x - x.^3));
    Yf = Yf .* H_bw;
    simResult.H_bw = H_bw;
    simResult.Yf = Yf;

    % Apply filter to signals
    simResult.output = real(ifft(Yf));
end

% function simResult = ApplyReceiverBW(simResult, simSettings)
% % GPU-accelerated 4th-order Butterworth low-pass filter
% 
%     % Backup previous output
%     simResult.outputPrevious = simResult.output;
% 
%     % Extract parameters
%     output = simResult.output;  % Assumed to be gpuArray
%     ts = simSettings.ts;
%     fc = simResult.receiverBW;
% 
%     % Design Butterworth filter (on CPU)
%     nyquist_freq = 1 / (2 * ts);
%     Wn = fc / nyquist_freq;
%     [b, a] = butter(4, Wn);
% 
%     % Move coefficients to GPU
%     b_gpu = gpuArray(b);
%     a_gpu = gpuArray(a);
% 
%     % Apply filter on GPU
%     simResult.output = filter(b_gpu, a_gpu, output);
% end
