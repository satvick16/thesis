function DisplayCTLE(simResult, simSettings)
% Displays the analog magnitude frequency response of the CTLE filter.
% X-axis is log-scaled and normalized to fb (bitrate), matching published figure style.
%
% INPUTS:
%   simResult.eq.ctle.g1, g2   - CTLE zero gains (dB)
%   simResult.eq.ctle.fp1-3    - CTLE pole frequencies (Hz)
%   simSettings.dataRate       - Data rate in Hz (used for normalization)

    % === Frequency vector ===
    f = simResult.fftFreqs(:);
    fb = simSettings.dataRate;

    % === Transfer functions ===
    H_ctle = simResult.H_ctle;
    H_ch = simResult.H_ch;

    % === Plot ===
    figure;
    h1 = semilogx(f / fb, 20*log10(abs(H_ctle)), 'LineWidth', 2, 'Color', 'Red');
    hold on;
    h2 = semilogx(f / fb, 20*log10(abs(H_ctle .* H_ch)), 'LineWidth', 2, 'Color', 'Blue');
    hold on;
    legend([h1, h2], {'CTLE', 'CTLE × Channel'}, 'Interpreter', 'latex');
    grid on;

    % === Mark poles ('x') and zeros ('o') on curve ===
    % Compute zero frequencies
    g1 = simResult.eq.ctle.g1;
    g2 = simResult.eq.ctle.g2;
    fp1 = simResult.eq.ctle.fp1;
    fp2 = simResult.eq.ctle.fp2;
    fp3 = simResult.eq.ctle.fp3;
    fz1 = 10^(g1/20) * fp1;
    fz2 = 10^(g2/20) * fp2;

    gain = 10^((g1 + g2)/20);
    wz1 = 2 * pi * fz1;
    wz2 = 2 * pi * fz2;
    wp1 = 2 * pi * fp1;
    wp2 = 2 * pi * fp2;
    wp3 = 2 * pi * fp3;

    % Interpolate magnitude at marker frequencies
    H_z1 = 20 * log10(abs( ...
        gain * (1 + 1j*2*pi*fz1/wz1) * (1 + 1j*2*pi*fz1/wz2) / ...
        ((1 + 1j*2*pi*fz1/wp1) * (1 + 1j*2*pi*fz1/wp2) * (1 + 1j*2*pi*fz1/wp3)) ));
    H_z2 = 20 * log10(abs( ...
        gain * (1 + 1j*2*pi*fz2/wz1) * (1 + 1j*2*pi*fz2/wz2) / ...
        ((1 + 1j*2*pi*fz2/wp1) * (1 + 1j*2*pi*fz2/wp2) * (1 + 1j*2*pi*fz2/wp3)) ));
    
    H_p1 = 20 * log10(abs( ...
        gain * (1 + 1j*2*pi*fp1/wz1) * (1 + 1j*2*pi*fp1/wz2) / ...
        ((1 + 1j*2*pi*fp1/wp1) * (1 + 1j*2*pi*fp1/wp2) * (1 + 1j*2*pi*fp1/wp3)) ));
    H_p2 = 20 * log10(abs( ...
        gain * (1 + 1j*2*pi*fp2/wz1) * (1 + 1j*2*pi*fp2/wz2) / ...
        ((1 + 1j*2*pi*fp2/wp1) * (1 + 1j*2*pi*fp2/wp2) * (1 + 1j*2*pi*fp2/wp3)) ));
    H_p3 = 20 * log10(abs( ...
        gain * (1 + 1j*2*pi*fp3/wz1) * (1 + 1j*2*pi*fp3/wz2) / ...
        ((1 + 1j*2*pi*fp3/wp1) * (1 + 1j*2*pi*fp3/wp2) * (1 + 1j*2*pi*fp3/wp3)) ));

    % Plot markers at true response height
    h_z1 = semilogx(fz1 / fb, H_z1, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5); % Zero 1
    h_z1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    h_z2 = semilogx(fz2 / fb, H_z2, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5); % Zero 2
    h_z2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    h_p1 = semilogx(fp1 / fb, H_p1, 'kx', 'MarkerSize', 8, 'LineWidth', 1.5); % Pole 1
    h_p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    h_p2 = semilogx(fp2 / fb, H_p2, 'kx', 'MarkerSize', 8, 'LineWidth', 1.5); % Pole 2
    h_p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    h_p3 = semilogx(fp3 / fb, H_p3, 'kx', 'MarkerSize', 8, 'LineWidth', 1.5); % Pole 3
    h_p3.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Final plot settings
    xlim([1e-3, 10]);
    ylim([-20, 0]);

    xlabel('Frequency (Normalized to $f_{b}$)', 'Interpreter', 'latex');
    ylabel('Magnitude (dB)', 'Interpreter', 'latex');
    title('CTLE Frequency Response', 'Interpreter', 'latex');

    set(gca, 'XTick', [1e-3 1e-2 1e-1 1 10]);
end
