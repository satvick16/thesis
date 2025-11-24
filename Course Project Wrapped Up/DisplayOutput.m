function DisplayOutput(simResult, simSettings, plotTitle)
% Displays a stylized plot of the output signal and cursor taps.
%
% INPUTS:
%   simResult.output           - Channel response signal
%   simResult.outputTime       - Time vector (in seconds)
%   simResult.cursorValues     - Tap amplitudes (in volts)
%   simResult.cursorTimes      - Tap locations (in seconds)
%   simSettings.ts             - Sampling interval (in seconds)
%   simSettings.samplesPerSymb - Samples per symbol (UI)
%   simSettings.plotPre        - Number of precursors to label
%   simSettings.plotPost       - Number of postcursors to label
%   plotTitle                  - Title of the figure

    % Extract and scale data
    output = simResult.output;
    ts = simSettings.ts;
    samplesPerSymb = simSettings.samplesPerSymb;

    outputTime = simResult.outputTime / (ts * samplesPerSymb);              % Time in UI
    cursorTimes = simResult.cursorTimes * 1e-9 / (ts * samplesPerSymb);     % Time in UI
    cursorValues = simResult.cursorValues * 1e3;                             % Amplitude in mV

    plotPre = simSettings.plotPre;
    plotPost = simSettings.plotPost;

    % Plot signal
    figure;
    plot(outputTime, output * 1e3, 'r', 'LineWidth', 2);  % Convert output to mV
    hold on;

    % Stem-like tap markers (hollow red circles)
    plot(cursorTimes, cursorValues, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);

    % Find and annotate main cursor (h_0)
    [~, mainIdx] = max(cursorValues);
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

    % Axes and styling
    xlabel('Time (UI)', 'Interpreter', 'latex');
    ylabel('Amplitude (mV)', 'Interpreter', 'latex');
    title(plotTitle, 'Interpreter', 'latex');
    grid on;

    % Set axis limits robustly
    leftBound = max(1, mainIdx - plotPre);
    rightBound = min(length(cursorTimes), mainIdx + plotPost);
    xlim([cursorTimes(leftBound) - 1, cursorTimes(rightBound) + 1]);

    ylim([-0.3 * max(abs(cursorValues)), 1.2 * max(cursorValues)]);

    hold off;
end
