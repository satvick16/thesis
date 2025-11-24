function simResult = ApplyPackage(simResult, simSettings)
% Applies package effect to the input pulse.
% INPUT:
%   simResult           - Struct containing input signal, fft frequencies
%   simSettings         - Struct with simulation settings
%
% OUTPUT:
%   simResult.input     - Filtered signal in time domain

    Yf = simResult.Yf;
    fftFreqs = simResult.fftFreqs(:);
    Cd = simSettings.Cd;
    Ls = simSettings.Ls;
    Rs = simSettings.Rs;
    Cb = simSettings.Cb;

    s = 1j*2*pi*fftFreqs;
    
    % Tranmitter package transfer function
    H_tx_package = 1 ./ (1 + (s.^2) * ((Ls+Rs) * Cb));

    % Receiver package transfer function 1/(1+s^2*(Ls*Cd))
    H_rx_package = 1 ./ (1 + (s.^2) * ((Ls+Rs) * Cd));
    
    Yf = Yf .* H_tx_package .* H_rx_package;
    simResult.H_tx_package = H_tx_package;
    simResult.H_rx_package = H_rx_package;
    simResult.Yf = Yf;
    simResult.output = real(ifft(Yf));
end