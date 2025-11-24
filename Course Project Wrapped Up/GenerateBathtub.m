function GenerateBathtub(simSettings, simResult)
% GenerateBathtub - Plot SER bathtub curves for each eye using 3D PDF
%
% Inputs:
%   - simSettings:   Struct with pdf3D, PAM level, pdf range, dx, etc.
%   - simResult:     Struct with .output used to extract mainCursor

    % === Axes ===
    samplesPerSymb = simSettings.samplesPerSymb;
    X = linspace(-0.5, 0.5, samplesPerSymb);        % Time axis (UI)
    Y = simSettings.pdfRange * 1e3;                 % Amplitude axis (mV)
    dx = simSettings.pdfdx * 1e3;                   % Bin width in mV
    pam = simSettings.pam;
    pdf3D = simResult.pdf3D;

    % === Compute true voltage levels from output ===
    mainCursor = max(simResult.output);
    levels = linspace(-mainCursor, mainCursor, pam) * 1e3;  % mV
    decisionBoundaries = 0.5 * (levels(1:end-1) + levels(2:end));  % mV

    % === Compute SER for each eye (region between boundaries) ===
    nEyes = pam - 1;
    SER = zeros(nEyes, samplesPerSymb);

    % Extend decision boundaries to full range
    eyeBoundaries = [-Inf, decisionBoundaries, Inf];
    
    for t = 1:samplesPerSymb
        pdf_slice = pdf3D(t, :) / (sum(pdf3D(t, :)) * dx);  % Normalize to unit mass
        cdf_slice = cumsum(pdf_slice) * dx;
    
        for eye = 1:nEyes
            lowB = eyeBoundaries(eye);
            highB = eyeBoundaries(eye+1);
    
            % Get bin indices in voltage axis
            lowIdx = find(Y >= lowB, 1, 'first');
            highIdx = find(Y >= highB, 1, 'first');
    
            % Guard against boundaries that fall outside Y range
            if isempty(lowIdx), lowIdx = 1; end
            if isempty(highIdx), highIdx = length(Y); end
    
            % SER = probability outside [lowB, highB]
            SER(eye, t) = 1 - (cdf_slice(highIdx) - cdf_slice(lowIdx));
        end
    end

    % === Plot Bathtub Curves ===
    colors = lines(nEyes);
    figure;
    for i = 1:nEyes
        subplot(nEyes, 1, i);
        semilogy(X, SER(i,:), 'LineWidth', 1.5, 'Color', colors(i,:));
        hold on;

        % Find opening at target SER (e.g. simSettings.targetSER)
        targetSER = simSettings.targetSER;
        idx = find(SER(i,:) < targetSER);
        if ~isempty(idx)
            left = X(idx(1));
            right = X(idx(end));
            opening = right - left;

            % Annotate opening
            plot([left, right], [targetSER, targetSER], 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
            text(mean([left, right]), targetSER * 2, sprintf('%.3f UI', opening), ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'Interpreter', 'latex');
        end

        ylabel('SER');
        ylim([1e-6, 1]);
        xlim([-0.5, 0.5]);
        grid on;
        legend(sprintf('Eye %d', i), 'Location', 'northeast', 'Interpreter', 'latex');
        if i == nEyes
            xlabel('Eye Opening (UI)');
        end
    end
end
