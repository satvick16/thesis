function DisplayOutputLegacy(simResult, simSettings, plotTitle)
% Displays a combined plot of the output signal and cursor values.
%
% INPUTS:
%   simResult.output              - Channel response signal
%   simResult.outputTime          - Time vector corresponding to the output
%   simResult.cursorValues        - Array containing the corresponding tap amplitudes
%   simResult.cursorTimes         - Array containing the extracted cursor times (in ns)
%   plotPre                       - Number of precursor cursors to plot
%   plotPost                      - Number of postcursor cursors to plot
%   plotTitle                     - Title of the figure
    
    output = simResult.output;
    outputTime = simResult.outputTime;
    cursorValues = simResult.cursorValues;
    cursorTimes = simResult.cursorTimes;
    plotPre = simSettings.plotPre;
    plotPost = simSettings.plotPost;
    
    % Stem plot for cursor values
    figure;
    plot(outputTime * 1e9, output, 'LineWidth', 2);
    hold on;
    
    % Plot Stem Data
    stem(cursorTimes, cursorValues, 'filled', 'LineWidth', 2);

    % Plot details
    title(plotTitle);
    xlabel('Time (ns)');
    ylabel('Amplitude');
    grid on;

    % Annotate Key Points
    [mainVal, mainIdx]= max(cursorValues);
    % Highlight Main Cursor
    stem(cursorTimes(mainIdx), mainVal, 'filled', 'LineWidth', 2, 'Color', 'k'); 
    text(cursorTimes(mainIdx), mainVal, ' h_0', 'VerticalAlignment', 'bottom');

    % Annotate Precursors
    for i = 1:plotPre
        text(cursorTimes(mainIdx-i), cursorValues(mainIdx-i), [' h_{-', num2str(i),'}'], 'VerticalAlignment', 'bottom');
    end

    % Annotate Postcursors
    for i = 1:plotPost
        text(cursorTimes(mainIdx+i), cursorValues(mainIdx+i), [' h_{', num2str(i),'}'], 'VerticalAlignment', 'bottom');
    end
    
    axis([cursorTimes(mainIdx-plotPre), cursorTimes(mainIdx+plotPost), -mainVal*0.5, mainVal*1.3]);
    hold off;
end
