function simResult = ApplyDriverFilter(simResult, simSettings)
% Applies the (Gaussian) driver filter to the input pulse.
% INPUT:
%   simResult           - Struct containing input signal in frequency domain, fft frequencies
%   simSettings         - Struct with simulation settings (dataRate, numSamples, Tt)
%
% OUTPUT:
%   simResult.output     - Filtered signal in time domain

    Yf = simResult.Yf;
    fftFreqs = simResult.fftFreqs;
    Tt = simSettings.Tt;
    dataRate = simSettings.dataRate;
    
    a = 0.8 / (Tt / dataRate);      % a = 0.8 / Tt (20%-80% rise/fall time in second)
    H_df = exp(-(pi * fftFreqs / a).^2);
    
    Yf = Yf .* H_df(:);
    
    simResult.H_df = H_df;
    simResult.Yf = Yf;
    simResult.output = real(ifft(Yf));
end