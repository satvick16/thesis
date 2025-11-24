function DisplayOutputComparison(simResult, simSettings, plotTitle, label1, label2)
% Plots the channel output and a previous signal (e.g., FFE output) in a
% stylized format similar to published figures.
%
% INPUTS:
%   simResult.output             - Current signal
%   simResult.outputPrevious     - Previous signal
%   simResult.outputTime         - Time vector (in seconds)
%   simResult.cursorValues       - Tap amplitudes (in volts)
%   simResult.cursorTimes        - Tap locations (in seconds)
%   simSettings.ts               - Sampling interval (in seconds)
%   simSettings.samplesPerSymb   - Samples per symbol (UI)
%   simSettings.plotPre          - Number of precursors to label
%   simSettings.plotPost         - Number of postcursors to label
%   plotTitle                    - Title of the figure
%   label1                       - Label 1 of legend (previous signal)
%   label2                       - Label 2 of legend (current signal)

    % Convert units
    ts = simSettings.ts;
    sps = simSettings.samplesPerSymb;

    timeUI = simResult.outputTime / (ts * sps);
    outputRed = simResult.output * 1e3;                 % Current signal (mV)
    outputBlack = simResult.outputPrevious * 1e3;       % Previous signal (mV)

    cursorTimes = simResult.cursorTimes * 1e-9 / (ts * sps); % Tap locations in UI
    cursorValues = simResult.cursorValues * 1e3;             % Tap amplitudes in mV

    plotPre = simSettings.plotPre;
    plotPost = simSettings.plotPost;

    % Main cursor index
    [~, mainIdx] = max(cursorValues);

    % Start plotting
    figure;

    % Plot previous signal (black line)
    plot(timeUI, outputBlack, 'k', 'LineWidth', 2); 
    hold on;

    % Plot current signal (red line)
    plot(timeUI, outputRed, 'r', 'LineWidth', 2); 

    % Hollow red circles for cursor values
    plot(cursorTimes, cursorValues, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);

    % Highlight main cursor with filled black circle
    % plot(cursorTimes(mainIdx), cursorValues(mainIdx), 'ko', 'MarkerSize', 8, 'LineWidth', 2);

    % Annotate h_0
    text(cursorTimes(mainIdx), cursorValues(mainIdx), '$h_0$', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

    % Annotate precursors
    for i = 1:plotPre
        idx = mainIdx - i;
        if idx >= 1
            text(cursorTimes(idx), cursorValues(idx), ...
                ['$h_{-', num2str(i), '}$'], ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
        end
    end

    % Annotate postcursors
    for i = 1:plotPost
        idx = mainIdx + i;
        if idx <= length(cursorValues)
            text(cursorTimes(idx), cursorValues(idx), ...
                ['$h_{', num2str(i), '}$'], ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
        end
    end

    % Labels, grid, and axis formatting
    xlabel('Time (UI)', 'Interpreter', 'latex');
    ylabel('Amplitude (mV)', 'Interpreter', 'latex');
    title(plotTitle, 'Interpreter', 'latex');
    legend(label1, label2, 'Location', 'northeast', 'Interpreter', 'latex');
    grid on;

    % Set axis bounds
    leftBound = max(1, mainIdx - plotPre);
    rightBound = min(length(cursorTimes), mainIdx + plotPost);
    yBound = max([max(outputBlack), max(abs(cursorValues))]);
    xlim([cursorTimes(leftBound) - 1, cursorTimes(rightBound) + 1]);
    ylim([-0.3 * yBound, 1.2 * yBound]);

    hold off;
end