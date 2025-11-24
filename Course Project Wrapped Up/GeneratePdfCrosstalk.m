function simResult = GeneratePdfCrosstalk(simResult, simSettings, pltSettings, simResources)
% Generates Crosstalk PDF by computing the output of each aggressor channel
% using its transfer function and applying the RX chain post-optimization.
%
% INPUTS:
%   simResult     - Simulation result
%   simSettings   - Simuation settings
%   pltSettings   - Plot settings updated from simSettings and simResult
%   simResources  - Struct with crosstalk rational models (rS11, rS12, rS21, rS22)
%
% OUTPUT:
%   pdfCrosstalk  - Final Crosstalk PDF from all aggressors
    
    fftFreqs = simResult.fftFreqs;
    
    N = simSettings.numSamples;
    pam = simSettings.pam;
    samplesPerSymb = simSettings.samplesPerSymb;

    pdfRange = pltSettings.pdfRange;
    dx = pltSettings.pdfdx;

    zs = simSettings.zs;
    zl = simSettings.zl;
    z0 = simSettings.z0;  % Differential impedance

    Gamma1 = (2 * zs - 2 * z0) / (2 * zs + 2 * z0);
    Gamma2 = (2 * zl - 2 * z0) / (2 * zl + 2 * z0);

    % === Apply Transceiver Signal Processing Chain ===
    tmpResult = simResult;
    Xf = fft(tmpResult.input, N);

    % === Initialize with Dirac delta ===
    pdfCrosstalk = zeros(1, length(pdfRange));
    [~, zeroIdx] = min(abs(pdfRange));
    pdfCrosstalk(zeroIdx) = 1 / dx;

    % === Loop over each crosstalk aggressor ===
    for i = 1:numel(simResources.crosstalk)
        % --- Extract S-parameter rational functions ---
        rS11 = simResources.crosstalk(i).rS11;
        rS12 = simResources.crosstalk(i).rS12;
        rS21 = simResources.crosstalk(i).rS21;
        rS22 = simResources.crosstalk(i).rS22;

        % --- Evaluate frequency-domain transfer function Hf ---
        S11 = squeeze(freqresp(rS11, fftFreqs));
        S12 = squeeze(freqresp(rS12, fftFreqs));
        S21 = squeeze(freqresp(rS21, fftFreqs));
        S22 = squeeze(freqresp(rS22, fftFreqs));

        DeltaS = S11 .* S22 - S12 .* S21;
        numerator = S21 .* (1 - Gamma1) .* (1 + Gamma2);
        denominator = 1 - S11 .* Gamma1 - S22 .* Gamma2 + Gamma1 .* Gamma2 .* DeltaS;
        H_xt = numerator ./ denominator;

        % --- Get time-domain crosstalk signal ---
        Yf = H_xt(:) .* Xf;
        Yf = Yf .* tmpResult.H_df;
        Yf = Yf .* tmpResult.H_tx_package;
        Yf = Yf .* tmpResult.H_rx_package;
        Yf = Yf .* tmpResult.H_bw;
        Yf = Yf .* tmpResult.H_ctle;
        Yf = Yf .* tmpResult.H_vga;
        tmpResult.output = real(ifft(Yf));
        tmpResult = ApplyFFE(tmpResult, tmpResult.eq.tapsTXFFE, simSettings);
        tmpResult = ApplyFFE(tmpResult, tmpResult.eq.tapsRXFFE, simSettings);

        % --- Construct PDF via oversampled Gaussian kernel ---
        partialPDFs = zeros(samplesPerSymb, length(pdfRange));
        for phase = 0:samplesPerSymb-1
            downsampled = tmpResult.output(phase+1:samplesPerSymb:end);
            sigma = sqrt(sum(downsampled.^2) * 1/3 * (pam+1)/(pam-1));
            partialPDFs(phase+1, :) = normpdf(pdfRange, 0, sigma);
        end

        % --- Average over all phases and convolve into cumulative PDF ---
        pdfAggressor = mean(partialPDFs, 1);
        pdfCrosstalk = conv(pdfCrosstalk, pdfAggressor, 'same') * dx;
    end

    simResult.pdfCrosstalk = pdfCrosstalk;
end

% function pdfCrosstalk = GeneratePdfCrosstalk(simResult, pltSettings, simResources)
% % Generates Crosstalk PDF by downsampling each aggressor signal, averaging
% % phase-shifted Gaussian PDFs, and convolving them using Dirac-delta initialization.
% %
% % INPUTS:
% %   simResult.outputCrosstalk - Cell array of crosstalk signals
% %   pltSettings.samplesPerSymb - Oversampling factor
% %   pltSettings.pdfRange       - X-axis vector for PDF
% %   pltSettings.pdfdx          - Bin width
% %
% % OUTPUT:
% %   pdfCrosstalk              - Final Crosstalk PDF from all aggressors
% 
%     outputCrosstalk = simResult.outputCrosstalk;
%     samplesPerSymb = pltSettings.samplesPerSymb;
%     pdfRange = pltSettings.pdfRange;
%     dx = pltSettings.pdfdx;
% 
%     numAggressors = length(outputCrosstalk);
% 
%     % === Initialize with Dirac delta at x=0 ===
%     pdfCrosstalk = zeros(1, length(pdfRange));
%     [~, zeroIdx] = min(abs(pdfRange));  % Find index closest to 0
%     pdfCrosstalk(zeroIdx) = 1 / dx;     % Discrete delta function
% 
%     for i = 1:numAggressors
%         xt = outputCrosstalk{i};
%         partialPDFs = zeros(samplesPerSymb, length(pdfRange));
% 
%         for phase = 0:samplesPerSymb-1
%             % Downsample with phase offset
%             downsampled = xt(phase+1:samplesPerSymb:end);
% 
%             % Compute sigma as geometric average
%             % sigma = rms(downsampled);
%             sigma = sqrt(sum(downsampled.^2));
% 
%             % Gaussian PDF centered at 0
%             partialPDFs(phase+1, :) = normpdf(pdfRange, 0, sigma);
%         end
% 
%         % Average PDFs from all phase offsets
%         pdfAggressor = mean(partialPDFs, 1);
% 
%         % Convolve with cumulative result
%         pdfCrosstalk = conv(pdfCrosstalk, pdfAggressor, 'same') * dx;
%     end
% end
