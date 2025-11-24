function simResult = ApplyCTLE(simResult)
% Applies a Continuous-Time Linear Equalizer (CTLE) in frequency domain.
% Zeros are placed based on gain shaping relative to pole locations.
%
% INPUTS:
%   simResult.Yf                 - Input signal in frequency domain
%   simResult.eq.ctle.g1, g2     - CTLE zero gains (dB)
%   simResult.eq.ctle.fp1-3      - CTLE pole frequencies (Hz)
%
% OUTPUT:
%   simResult.H_ctle             - CTLE transfer function
%   simResult.Yf                 - CTLE-equalized signal in frequency domain
%   simResult.output             - CTLE-equalized signal in time domain

    % === Save current output before CTLE ===
    simResult.outputPrevious = simResult.output;

    % === Inputs ===
    fftFreqs = simResult.fftFreqs;
    Yf = simResult.Yf;
    g1 = simResult.eq.ctle.g1;
    g2 = simResult.eq.ctle.g2;
    fp1 = simResult.eq.ctle.fp1;
    fp2 = simResult.eq.ctle.fp2;
    fp3 = simResult.eq.ctle.fp3;

    % === Sanity Check ===
    if any([fp1, fp2, fp3] <= 0) || any(isnan([fp1, fp2, fp3])) ...
       || any(isnan([g1, g2])) || any(abs([g1, g2]) > 80)
        warning('CTLE skipped due to invalid parameter values.');
        return;
    end

    w = 2 * pi * fftFreqs;

    % === Compute CTLE analog transfer function ===
    % Place zeros using gain shaping relative to pole locations
    ctleGain = 10^((g1 + g2)/20);
    fz1 = fp1 * 10^(g1/20);
    fz2 = fp2 * 10^(g2/20);
    wz1 = 2 * pi * fz1;
    wz2 = 2 * pi * fz2;

    wp1 = 2 * pi * fp1;
    wp2 = 2 * pi * fp2;
    wp3 = 2 * pi * fp3;

    % Transfer function H(jw)
    H_ctle = ctleGain .* (1 + 1j*w/wz1) .* (1 + 1j*w/wz2) ./ ...
                         ((1 + 1j*w/wp1) .* (1 + 1j*w/wp2) .* (1 + 1j*w/wp3));
    simResult.H_ctle = H_ctle;

    % === Apply CTLE in Frequency Domain ===
    Yf = Yf .* H_ctle;
    simResult.Yf = Yf;
    simResult.output = real(ifft(Yf));
end

% function simResult = ApplyCTLE(simResult, simSettings)
% % GPU-accelerated CTLE application in the frequency domain
% 
%     % === Save current output before CTLE ===
%     simResult.outputPrevious = simResult.output;
% 
%     % === Inputs ===
%     output = simResult.output(:);  % Should already be gpuArray
%     g1 = simResult.eq.ctle.g1;
%     g2 = simResult.eq.ctle.g2;
%     fp1 = simResult.eq.ctle.fp1;
%     fp2 = simResult.eq.ctle.fp2;
%     fp3 = simResult.eq.ctle.fp3;
%     ts  = simSettings.ts;
%     numSamples = length(output);
%     N = 2^nextpow2(numSamples);
% 
%     % === Sanity Check ===
%     if any([fp1, fp2, fp3] <= 0) || any(isnan([fp1, fp2, fp3])) ...
%        || any(isnan([g1, g2])) || any(abs([g1, g2]) > 80)
%         warning('CTLE skipped due to invalid parameter values.');
%         return;
%     end
% 
%     % === Frequency Vector (on GPU) ===
%     f = gpuArray((0:N-1)' / (N * ts));  % Hz
%     w = 2 * pi * f;
% 
%     % === Compute CTLE transfer function on GPU ===
%     fz1 = fp1 * 10^(g1/20);
%     fz2 = fp2 * 10^(g2/20);
%     wz1 = 2 * pi * fz1;
%     wz2 = 2 * pi * fz2;
%     wp1 = 2 * pi * fp1;
%     wp2 = 2 * pi * fp2;
%     wp3 = 2 * pi * fp3;
% 
%     H = (1 + 1j*w/wz1) .* (1 + 1j*w/wz2) ./ ...
%         ((1 + 1j*w/wp1) .* (1 + 1j*w/wp2) .* (1 + 1j*w/wp3));
% 
%     % === Apply CTLE on GPU ===
%     OUTPUT = fft(output, N);
%     Y = OUTPUT .* H;
%     y = real(ifft(Y));
% 
%     % === Assign trimmed result back to simResult ===
%     simResult.output = y(1:numSamples);
% end