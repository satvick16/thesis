function SNR = EvalSNR(simResult,simSettings)
% EVAL_SNR Computes SNR from cursor values (post RX FFE, pre DFE)
%
% INPUT:
%   simResult.cursorValues      - Vector of sampled cursor amplitudes (precursors, main cursor, postcursors)
%                                 e.g., [h_{-3}, h_{-2}, h_{-1}, h_0, h_1, h_2, ...]
%
% OUTPUT:
%   SNR          - Signal-to-noise ratio computed as h0^2 / sum(h_i^2 for i ≠ 0)
    
    cursorValues = simResult.cursorValues;
    pam = simSettings.pam;

    % Identify index of main cursor (assumed to be the peak)
    [main_cursor, ~] = max(abs(cursorValues));

    % Calculate correction factor
    pamVarianceNormFactor = 1/3 * (pam + 1) / (pam - 1);

    % Signal power = main cursor squared
    signal_power = main_cursor^2 * pamVarianceNormFactor;

    % Residual power = sum of all other cursor energies
    residual_power = sum(cursorValues.^2) * pamVarianceNormFactor - signal_power;
    
    % Noise power = sum of tx and rx noise power
    [TXNoiseVar,RXNoiseVar] = ComputeNoiseVariance(simResult,simSettings);

    % Compute SNR
    SNR = 10 * log10(signal_power / (residual_power + TXNoiseVar + RXNoiseVar));
end


function [TXNoiseVar, RXNoiseVar] = ComputeNoiseVariance(simResult, simSettings)
% Computes the Variance of TX Noise and RX Noise.
%
% TX noise is shaped through: TX FFE → Channel → RX BW → CTLE → VGA → RX FFE
% RX noise is shaped through: RX BW → CTLE → VGA → RX FFE

    % === Parameters ===
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
    TXNoiseVar = txPSD * G_tx;
    RXNoiseVar = rxPSD * G_rx;
end
