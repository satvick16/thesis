function simResult = GeneratePdfJitter(simResult, simSettings, pltSettings)
% Generates jitter PDFs for each PAM transition level (outer to inner).
% Returns a cell array: pdfJitter{1}, pdfJitter{2}, ..., pdfJitter{pam/2}

    % === Extract settings ===
    pdfRange         = pltSettings.pdfRange;
    output           = simResult.output;
    jitterSampling   = simSettings.jitterSampling;
    UI               = 1 / simSettings.dataRate;
    ts               = simSettings.ts;
    diracSeparation  = simSettings.diracSeparation;
    samplesPerSymb   = simSettings.samplesPerSymb;
    trackPre         = simSettings.trackPre;
    trackPost        = simSettings.trackPost;
    pam              = simSettings.pam;

    % === Main tap ===
    [~, mainIdx] = max(output);

    % === Cursor indices ===
    totalTaps = trackPre + 1 + trackPost;
    center = trackPre + 1;
    cursorIndices = nan(1, totalTaps);
    cursorIndices(center) = mainIdx;

    for k = 1:trackPre
        idx = mainIdx - k * samplesPerSymb;
        if idx >= 1
            cursorIndices(center - k) = idx;
        end
    end
    for k = 1:trackPost
        idx = mainIdx + k * samplesPerSymb;
        if idx <= length(output)
            cursorIndices(center + k) = idx;
        end
    end
    cursorIndices = cursorIndices(~isnan(cursorIndices));

    % === Slopes ===
    d_output_dt = gradient(output, ts);
    baseSlope = sqrt(sum(abs(d_output_dt(cursorIndices)).^2));  % full-scale transition slope

    % === Level scaling factors for each unique transition (outer to inner) ===
    halfLevels = ceil(pam / 2);
    scaleVec = (pam - 1 - 2 * (0:halfLevels-1)) / (pam - 1);  % e.g., [1, 1/3] for PAM-4

    % === Precompute ===
    dx = mean(diff(pdfRange));
    pdfJitter = cell(1, halfLevels);

    for i = 1:halfLevels
        levelScale = scaleVec(i);

        % Voltage std for this transition level
        sigma_v = jitterSampling * UI * baseSlope * levelScale;

        % Gaussian PDF
        if sigma_v == 0
            [~, idx] = min(abs(pdfRange));
            pdfGauss = zeros(size(pdfRange)); pdfGauss(idx) = 1 / dx;
        else
            pdfGauss = normpdf(pdfRange, 0, sigma_v);
        end
        

        % Dirac delta shift (in voltage domain)
        eff_slope = baseSlope * levelScale;
        diracShiftVolt = (diracSeparation * UI / 2) * eff_slope;

        [~, idx1] = min(abs(pdfRange - (-diracShiftVolt)));
        [~, idx2] = min(abs(pdfRange - ( diracShiftVolt)));

        pdfDirac1 = zeros(size(pdfRange)); pdfDirac1(idx1) = 1 / dx;
        pdfDirac2 = zeros(size(pdfRange)); pdfDirac2(idx2) = 1 / dx;

        % Convolve: Gaussian * Dirac1 * Dirac2
        tmp = conv(pdfGauss, pdfDirac1, 'same') * dx;
        tmp = conv(tmp, pdfDirac2, 'same') * dx;

        % Normalize
        pdfJitter{i} = tmp / (sum(tmp) * dx);
    end

    simResult.pdfJitter = pdfJitter;
end
