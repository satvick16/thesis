function simResult = GeneratePdfNoise(simResult, simSettings, pltSettings)
% Generates TX and RX noise PDFs based on analytical shaping through transfer functions.
%
% TX noise is shaped through: TX FFE → Channel → RX BW → CTLE → VGA → RX FFE
% RX noise is shaped through: RX BW → CTLE → VGA → RX FFE
%
% INPUTS:
%   simResult     : Contains equalization settings (FFE, CTLE, VGA, etc.)
%   simSettings   : Contains pdfRange, SNDR (dB), rxPSD (V^2/GHz), dataRate, ts
%   pltSettings   : Plot Settings
%
% OUTPUTS:
%   pdfNoiseTX    : PDF of TX noise at decision point
%   pdfNoiseRX    : PDF of RX noise at decision point

    % === Parameters ===
    pdfRange = pltSettings.pdfRange;

    pam = simSettings.pam;
    txSNDR = simSettings.txSNDR;            % in dB
    rxPSD = simSettings.rxPSD * 1e-9;       % Convert from V^2/GHz to V^2/Hz
    f = simResult.fftFreqs(:);
    df = f(2) - f(1);

    % === TX Noise Power at Transmitter ===
    P_signal = (1/3) * (pam + 1) / (pam - 1);
    P_noise_tx = P_signal / (10^(txSNDR / 10));
    txPSD = P_noise_tx / 56e9;

    % === Channel Transfer Function ===
    H_ch = simResult.H_ch;

    % === RX BW (4th-Order Butterworth) ===
    H_bw = simResult.H_bw;

    % === CTLE Transfer Function ===
    H_ctle = simResult.H_ctle;

    % === VGA Transfer Function ===
    H_vga = simResult.H_vga;

    % === FFE taps ===
    h_tx = simResult.eq.tapsTXFFE(:).';
    h_rx = simResult.eq.tapsRXFFE(:).';

    % === Discrete-time FFE responses at symbol rate, mapped onto analog f ===
    Fs = simSettings.dataRate;
    w = 2*pi*f/Fs;
    E_tx = exp(-1j * w * (0:(numel(h_tx)-1)));    % TX FFE DT response over f
    E_rx = exp(-1j * w * (0:(numel(h_rx)-1)));    % RX FFE DT response over f
    
    H_tx_ffe = E_tx * h_tx.';
    H_rx_ffe = E_rx * h_rx.';

    % === Composite TX Noise Transfer Function ===
    H_tx = H_ch(:) .* H_bw(:) .* H_ctle(:) .* H_vga(:);
    G_tx = sum(abs(H_tx_ffe).^2 .* abs(H_tx).^2 .* abs(H_rx_ffe).^2) * df;

    % === Composite RX Noise Transfer Function ===
    H_rx = H_bw(:) .* H_ctle(:) .* H_vga(:);
    G_rx = sum(abs(H_rx).^2 .* abs(H_rx_ffe).^2) * df;

    % === Final Noise Standard Deviations ===
    sigma_tx = sqrt(txPSD * G_tx);
    sigma_rx = sqrt(rxPSD * G_rx);

    % === Generate PDFs ===
    pdfNoiseTX = normpdf(pdfRange, 0, sigma_tx);
    pdfNoiseRX = normpdf(pdfRange, 0, sigma_rx);

    simResult.pdfNoiseTX = pdfNoiseTX;
    simResult.pdfNoiseRX = pdfNoiseRX;
end
