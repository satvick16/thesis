function simResult = ApplyVGA(simResult)
% Applies a Variable Gain Amplifier (VGA) with a first-order low-pass filter.
%
% INPUTS:
%   simResult.Yf                 - Input signal in frequency domain
%   simResult.eq.vga.gain        - VGA Gain (dB)
%   simResult.eq.vga.bw          - 3dB bandwidth (Hz)
%
% OUTPUT:
%   simResult.H_vga              - Transfer function
%   simResult.Yf                 - Equalized signal in frequency domain
%   simResult.output             - Equalized signal in time domain

    % Backup previous output
    simResult.outputPrevious = simResult.output;

    % Extract parameters
    fftFreqs = simResult.fftFreqs(:);
    Yf = simResult.Yf;
    vgaGain = simResult.eq.vga.gain;
    vgaBW = simResult.eq.vga.bw;
    wc = 2 * pi * vgaBW;
    w = 2 * pi * fftFreqs;

    % Transfer function: H_vga(s) = K * wc / (s + wc)
    H_vga = 10^(vgaGain / 20) * (wc ./ (1j * w + wc));
    simResult.H_vga = H_vga;

    % Apply filter to main output and noise
    Yf = Yf .* H_vga;
    simResult.Yf = Yf;
    simResult.output = real(ifft(Yf));
end

% function simResult = ApplyVGA(simResult, simSettings)
% % GPU-accelerated VGA model using a first-order low-pass filter
% 
%     % Backup previous output
%     simResult.outputPrevious = simResult.output;
% 
%     % Extract parameters
%     output = simResult.output;  % Assumed to be gpuArray
%     gain = 10^(simResult.eq.vga.gain / 20);
%     bw = simResult.eq.vga.bw;
%     ts = simSettings.ts;
% 
%     % Continuous-time transfer function H(s) = K / (s + wc)
%     wc = 2 * pi * bw;
%     num = [gain * wc];
%     den = [1, wc];
% 
%     % Discretize via bilinear transform (still on CPU)
%     [b, a] = bilinear(num, den, 1 / ts);
% 
%     % Move to GPU
%     b_gpu = gpuArray(b);
%     a_gpu = gpuArray(a);
% 
%     % Filter on GPU
%     simResult.output = filter(b_gpu, a_gpu, output);
% end
