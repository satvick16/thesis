function DisplayPDF3D(simResult, simSettings, pltSettings)
% Displays 3D PDFs

    samplesPerSymb = simSettings.samplesPerSymb;
    dT = 1 / samplesPerSymb;
    pdf3D = simResult.pdf3D;
    SER = simSettings.targetSER;

    pdfRange = pltSettings.pdfRange;
    dV = pltSettings.pdfdx;

    % === Downsample pdf3D for plotting ===
    X = linspace(-0.5, 0.5, samplesPerSymb);
    Y = pdfRange * 1e3;
    downsampleFactor = 100;
    pdf3D_ds = pdf3D(:, 1:downsampleFactor:end);  % Still rows = samplesPerSymb
    Y_ds = Y(1:downsampleFactor:end);             % Downsample Y axis accordingly
    [X_mesh, Y_mesh] = meshgrid(X, Y_ds);

    % --- Plot 3D Surface Plot (Default View) ---
    figure;
    surf(X_mesh, Y_mesh, pdf3D_ds', 'EdgeColor', 'none');
    title('3-D PDF Plot (Default View)', 'Interpreter', 'latex');
    xlabel('Time Points (UI)', 'Interpreter', 'latex'); 
    ylabel('Amplitude (mV)', 'Interpreter', 'latex');
    zlabel('Probability Density', 'Interpreter', 'latex');
    ylim([min(Y)*0.2,max(Y)*0.2]);
    set(gca, 'XTick', [-0.5 -0.25 0 0.25 0.5]);
    view(3);

    % --- Plot Top-Down View (Heatmap) ---
    figure;
    surf(X_mesh, Y_mesh, pdf3D_ds', 'EdgeColor', 'none');
    colormap("hot");
    title('3-D PDF Plot (Heatmap)', 'Interpreter', 'latex');
    xlabel('Time Points (UI)', 'Interpreter', 'latex');
    ylabel('Amplitude (mV)', 'Interpreter', 'latex');
    zlabel('Probability Density', 'Interpreter', 'latex');
    ylim([min(Y)*0.2,max(Y)*0.2]);
    set(gca, 'XTick', [-0.5 -0.25 0 0.25 0.5]);
    view(2);

    % === Plot the SER contour ===
    figure;
    hold on;
    contour(X, Y_ds, pdf3D_ds', 100);  % 100 levels

    % --- Create a masked copy of the 3D pdf for contour extraction ---
    mainCursor = max(simResult.cursorValues) * 1e3;  % In mV
    keepIdx = abs(Y_ds) <= mainCursor;
    pdf3D_masked = pdf3D_ds(:, keepIdx);
    Y_trunc = Y_ds(keepIdx);

    % --- Compute z-threshold on full pdf ---
    z_thr = FindSERContourThreshold(pdf3D, dT, dV, SER);

    % --- Plot SER contour in eye only
    contour(X, Y_trunc, pdf3D_masked', [z_thr, z_thr], 'LineWidth', 2, 'LineColor', 'r');
    colormap('parula')
    xlabel('Time (UI)', 'Interpreter', 'latex');
    ylabel('Amplitude (mV)', 'Interpreter', 'latex');
    ser_str = sprintf('%.1f \\times 10^{%d}', SER / 10^floor(log10(SER)), floor(log10(SER)));
    title(['3-D Eye Density Contour (Target SER = $' ser_str '$)'], 'Interpreter', 'latex');
    set(gca, 'XTick', [-0.5 -0.25 0 0.25 0.5]);
    yRange = [min(Y_ds)*0.2, max(Y_ds)*0.2];
    ylim(yRange);
    xlim([-0.5, 0.5]);

    % === Plot the Bathtub curves ===
    bt = GenerateBathtubCurves(pdf3D, X, Y, dT, dV, mainCursor, simSettings.pam, SER, 100);

    nEyes = simSettings.pam - 1;
    colors = lines(nEyes);
    
    % --- Horizontal Bathtub Curve ---
    figure;
    for i = 1:nEyes
        subplot(nEyes, 1, i);
        hold on;
    
        left  = bt.horizontal.left{i};
        right = bt.horizontal.right{i};
        ser   = bt.ser_values;
    
        semilogy(left, ser, '-', 'LineWidth', 2, 'Color', colors(i,:));   % left edge
        semilogy(right, ser, '-', 'LineWidth', 2, 'Color', colors(i,:));  % right edge
    
        % --- Annotate eye opening at target SER ---
        [~, idx] = min(abs(ser - SER));
        xLeft = left(idx);
        xRight = right(idx);
        y = ser(idx);
        opening = xRight - xLeft;
        plot([xLeft, xRight], [y, y], 'k-', 'LineWidth', 1);  % black horizontal line
    
        % Plot and annotate
        semilogy([xLeft, xRight], [ser(idx), ser(idx)], 'ko', 'MarkerFaceColor', 'w');
        text(mean([xLeft, xRight]), ser(idx)*1.2, sprintf('%.3f UI\n', opening), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'k', 'Interpreter', 'latex');
    
        set(gca, 'YScale', 'log');
        xlim([-0.5, 0.5]);
        ylim([SER / 100, 1]);
        grid on;
    
        ylabel('SER', 'Interpreter', 'latex');
        legend({sprintf('Eye %d', i)}, 'Location', 'southeast', 'Interpreter', 'latex');
    
        if i == nEyes
            xlabel('Eye Opening (UI)', 'Interpreter', 'latex');
        end
    end
    sgtitle('Horizontal Bathtub Curve', 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    % --- Vertical Bathtub Curve ---
    figure;
    hold on;
    hLegends = gobjects(nEyes, 1);  % Handles for legend
    
    for i = 1:nEyes
        top    = bt.vertical.top{i};
        bottom = bt.vertical.bottom{i};
        ser    = bt.ser_values;
    
        % Capture top edge line for legend
        hLegends(i) = semilogx(ser, top, '-', 'LineWidth', 2, 'Color', colors(i,:));
        semilogx(ser, bottom, '-', 'LineWidth', 2, 'Color', colors(i,:));
    
        % --- Annotation ---
        [~, idx] = min(abs(ser - SER));
        yTop = top(idx);
        yBottom = bottom(idx);
        x = ser(idx);
        opening = yTop - yBottom;
    
        plot([x, x], [yBottom, yTop], 'k-', 'LineWidth', 1);
        semilogx(ser(idx), [yBottom, yTop], 'ko-', 'MarkerFaceColor', 'w');
    
        % Text
        text(ser(idx)*1.1, mean([yTop, yBottom]), sprintf('%.1f mV', opening), ...
            'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', 'k', 'Interpreter', 'latex');
    end
    
    % Final plot setup
    set(gca, 'XScale', 'log');
    xlabel('SER', 'Interpreter', 'latex');
    ylabel('Eye Opening (mV)', 'Interpreter', 'latex');
    ylim([-mainCursor*1.2, mainCursor*1.2]);
    xlim([SER / 100, 1]);
    legend(hLegends, {'Top Eye', 'Middle Eye', 'Bottom Eye'}, 'Location', 'northwest', 'Interpreter', 'latex');
    grid on;
    title('Vertical Bathtub Curve', 'Interpreter', 'latex');
end

function z_threshold = FindSERContourThreshold(pdf3D, dx, dy, targetSER)
    % Flatten and sort all PDF values in ascending order
    z = sort(pdf3D(:), 'ascend');

    % Compute volume per bin
    dz = dx * dy;

    % Cumulative integral of increasing density
    cumulative_mass = cumsum(z) * dz;

    % Find first index where mass exceeds SER
    idx = find(cumulative_mass >= targetSER, 1, 'first');

    % Threshold is the z-value where the accumulated tail equals SER
    z_threshold = z(idx);
end

function bt = GenerateBathtubCurves(pdf3D, X, Y, dx, dy, mainCursor, pam, targetSER, numLevels)
    % Precompute constants
    serLevels = logspace(log10(targetSER / 100), log10(0.5), numLevels);
    numEyes = pam - 1;
    levels = linspace(-mainCursor, mainCursor, pam);

    % Preallocate outputs as flat arrays
    left  = nan(numEyes, numLevels);
    right = nan(numEyes, numLevels);
    top   = nan(numEyes, numLevels);
    bottom = nan(numEyes, numLevels);

    % === Parallel Loop ===
    parfor i = 1:numLevels
        thr = FindSERContourThreshold(pdf3D, dx, dy, serLevels(i));
        M = contourc(X, Y, pdf3D', [thr, thr]);

        % Temporary storage for each eye at this SER
        left_i  = nan(1, numEyes);
        right_i = nan(1, numEyes);
        top_i   = nan(1, numEyes);
        bottom_i = nan(1, numEyes);

        idx = 1;
        while idx < size(M,2)
            numPoints = M(2, idx);
            xPoints = M(1, idx+1 : idx+numPoints);
            yPoints = M(2, idx+1 : idx+numPoints);
            idx = idx + numPoints + 1;

            for eye = 1:numEyes
                yLow = levels(eye);
                yHigh = levels(eye+1);
                inBand = yPoints >= yLow & yPoints <= yHigh;

                if any(inBand)
                    xEye = xPoints(inBand);
                    yEye = yPoints(inBand);

                    left_i(eye)   = min([left_i(eye), xEye], [], 'omitnan');
                    right_i(eye)  = max([right_i(eye), xEye], [], 'omitnan');
                    bottom_i(eye) = min([bottom_i(eye), yEye], [], 'omitnan');
                    top_i(eye)    = max([top_i(eye), yEye], [], 'omitnan');
                end
            end
        end

        % Store back into output matrices
        left(:, i)   = left_i;
        right(:, i)  = right_i;
        top(:, i)    = top_i;
        bottom(:, i) = bottom_i;
    end

    % === Reshape results into struct-of-cells ===
    for eye = 1:numEyes
        bt.horizontal.left{eye}  = left(eye, :);
        bt.horizontal.right{eye} = right(eye, :);
        bt.vertical.top{eye}     = top(eye, :);
        bt.vertical.bottom{eye}  = bottom(eye, :);
    end
    bt.ser_values = serLevels;
end
