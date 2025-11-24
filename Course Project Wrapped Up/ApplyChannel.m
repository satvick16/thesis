function simResult = ApplyChannel(simResult, simSettings, simResources)
% Applies channel filtering using precomputed rational models stored in simResources.
% INPUT:
%   simResult           - Struct containing output signal
%   simSettings         - Struct with simulation settings (ts, zs, zl, z0, numSamples)
%   simResources        - Struct with pre-fitted rational S-parameter models
%
% OUTPUT:
%   simResult.Yf        - Output Signal in Frequency Domain
%   simResult.output    - Output Signal in Time Domain

    % Backup current output
    simResult.outputPrevious = simResult.output;

    Yf = simResult.Yf;
    zs = simSettings.zs;
    zl = simSettings.zl;
    fftFreqs = simResult.fftFreqs;

    % Retrieve rational S-parameter models
    rS11 = simResources.channel.rS11;
    rS12 = simResources.channel.rS12;
    rS21 = simResources.channel.rS21;
    rS22 = simResources.channel.rS22;

    % Evaluate frequency response
    S11 = squeeze(freqresp(rS11, fftFreqs));
    S12 = squeeze(freqresp(rS12, fftFreqs));
    S21 = squeeze(freqresp(rS21, fftFreqs));
    S22 = squeeze(freqresp(rS22, fftFreqs));

    % Compute reflection coefficients
    z0 = 2 * simSettings.z0;
    Gamma1 = (2 * zs - z0) / (2 * zs + z0);
    Gamma2 = (2 * zl - z0) / (2 * zl + z0);

    % Compute H_ch
    DeltaS = S11 .* S22 - S12 .* S21;
    numerator = S21 .* (1 - Gamma1) .* (1 + Gamma2);
    denominator = 1 - S11 .* Gamma1 - S22 .* Gamma2 + Gamma1 * Gamma2 * DeltaS;
    H_ch = numerator ./ denominator;
    simResult.H_ch = H_ch;

    % Apply channel to input signal
    Yf = H_ch(:) .* Yf;
    simResult.Yf = Yf;
    simResult.output = real(ifft(Yf));
end

% function simResult = ApplyChannel(simResult, simSettings, simResources)
% % GPU-accelerated ApplyChannel
% % Applies channel filtering using rational models and GPU-based FFT
% 
%     input = simResult.input;     % Already assumed to be on GPU (gpuArray)
%     ts = simSettings.ts;
%     zs = simSettings.zs;
%     zl = simSettings.zl;
%     numSamples = simSettings.numSamples;
% 
%     % Generate FFT frequency vector on GPU
%     fftFreqs = gpuArray((0:numSamples-1)') / (ts * numSamples);
% 
%     % Retrieve rational S-parameter models
%     rS11 = simResources.channel.rS11;
%     rS12 = simResources.channel.rS12;
%     rS21 = simResources.channel.rS21;
%     rS22 = simResources.channel.rS22;
% 
%     % Evaluate frequency response (still CPU-side since freqresp doesn't support gpuArray)
%     % But this is fast compared to FFT
%     S11 = gpuArray(squeeze(freqresp(rS11, fftFreqs)));
%     S12 = gpuArray(squeeze(freqresp(rS12, fftFreqs)));
%     S21 = gpuArray(squeeze(freqresp(rS21, fftFreqs)));
%     S22 = gpuArray(squeeze(freqresp(rS22, fftFreqs)));
% 
%     % Compute reflection coefficients
%     z0 = 2 * simSettings.z0;
%     Gamma1 = (2 * zs - z0) / (2 * zs + z0);
%     Gamma2 = (2 * zl - z0) / (2 * zl + z0);
% 
%     % Compute H(f)
%     DeltaS = S11 .* S22 - S12 .* S21;
%     numerator = S21 .* (1 - Gamma1) .* (1 + Gamma2);
%     denominator = 1 - S11 .* Gamma1 - S22 .* Gamma2 + Gamma1 * Gamma2 * DeltaS;
%     Hf = numerator ./ denominator;
% 
%     % Backup current output
%     simResult.outputPrevious = input;
% 
%     % Apply channel to input signal in frequency domain on GPU
%     Xf = fft(input(:), numSamples);
%     Yf = Hf(:) .* Xf;
%     simResult.output = real(ifft(Yf));
% end