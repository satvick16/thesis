function DisplayPDF2D(simResult, pltSettings)
% Displays 2D PDFs in a 2x4 layout where the total signal spans both rows of the rightmost column.

    pdfRange = pltSettings.pdfRange;
    amplitude = pdfRange * 1e3;

    % Extract PDFs
    pdfMainCursor = simResult.pdfMainCursor;
    pdfISI        = simResult.pdfISI;
    pdfCrosstalk  = simResult.pdfCrosstalk;
    pdfNoiseTX    = simResult.pdfNoiseTX;
    pdfNoiseRX    = simResult.pdfNoiseRX;
    pdfJitter     = simResult.pdfJitter;
    pdfTotal      = simResult.pdfTotal;

    % Create figure
    figure('Color', 'w', 'Position', [100 100 1400 600]);

    % Subplots (a)–(f)
    pdfList = {pdfMainCursor, pdfISI, pdfCrosstalk, pdfNoiseTX, pdfNoiseRX};
    titles  = {'Main Cursor', 'Residual ISI', 'Crosstalk', 'TX Noise', 'RX Noise'};
    labels  = {'(a)', '(b)', '(c)', '(d)', '(e)'};

    subplotIdx = [1 2 3 5 6];  % MATLAB subplot positions (2 rows × 4 columns)
    for i = 1:5
        subplot(2, 4, subplotIdx(i));
        plot(pdfList{i}, amplitude, 'LineWidth', 1.5);
        ylabel('Amplitude (mV)');
        title(sprintf('%s %s', labels{i}, titles{i}));
        xlim([0,max(pdfList{i})*1.2]);
        if i~=1
            ylim([min(amplitude)*0.05,max(amplitude)*0.05]);
        else
            ylim([min(amplitude)*0.15,max(amplitude)*0.15]);
        end
        set(gca, 'FontSize', 10, 'LineWidth', 1, 'YDir', 'normal', 'XColor', [1 1 1], 'GridColor', [0 0 0]);
        grid on;
    end

    % Subplot (f) Jitter — position (7)
    subplot(2,4,7);
    hold on;
    numLevels = numel(pdfJitter);
    colors = lines(numLevels);

    for i = 1:numLevels
        plot(pdfJitter{i}, amplitude, 'Color', colors(i,:), 'LineWidth', 1.5);
    end

    % Jitter legend
    switch numLevels
        case 1, labels = {'Outer Levels'};
        case 2, labels = {'Outer Levels', 'Inner Levels'};
        case 3, labels = {'Outer', 'Mid', 'Inner'};
        case 4, labels = {'Outer', 'Mid-1', 'Mid-2', 'Inner'};
        otherwise
            labels = arrayfun(@(i) sprintf('Level %d', i), 1:numLevels, 'UniformOutput', false);
    end
    legend(labels, 'FontSize', 9, 'Location', 'southwest', 'Interpreter', 'latex');

    ylabel('Amplitude (mV)', 'Interpreter', 'latex');
    title('(f) Jitter Noise', 'Interpreter', 'latex');
    xlim([0, max(cellfun(@max, pdfJitter))*1.2]);
    ylim([min(amplitude)*0.05,max(amplitude)*0.05]);
    set(gca, 'FontSize', 10, 'LineWidth', 1, 'YDir', 'normal', 'XColor', [1 1 1], 'GridColor', [0 0 0]);
    grid on;
    hold off;

    % Subplot (g) Total Signal — spanning two rows in column 4
    subplot(2,4,[4 8]);
    plot(pdfTotal, amplitude, 'r', 'LineWidth', 2);
    ylabel('Amplitude (mV)', 'Interpreter', 'latex');
    title('(g) Total Signal', 'Interpreter', 'latex');
    xlim([0,max(pdfTotal)*1.2]);
    ylim([min(amplitude)*0.15,max(amplitude)*0.15]);
    set(gca, 'FontSize', 10, 'LineWidth', 1, 'YDir', 'normal', 'XColor', [1 1 1], 'GridColor', [0 0 0]);
    grid on;

    % Global title
    sgtitle(sprintf('Probability Density Functions (PDFs) at Slicer Input\n'), ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    fprintf('Estimated SER: %.4e\n', simResult.SER);
end
