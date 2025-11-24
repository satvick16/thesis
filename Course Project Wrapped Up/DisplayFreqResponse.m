function DisplayFreqResponse(simSettings)
%   Plots frequency responses of the channel and crosstalk files
%
%   Inputs:
%       simSettings.channelFilename         - Main channel S-parameter file (.s4p)
%       simSettings.crosstalk_filenames     - Cell array of crosstalk S-parameter files (.s4p)
%       simSettings.zl                      - Load impedance in ohms
%       simSettings.zs                      - Source impedance in ohms
    
    channelFilename = simSettings.channelFilename;
    crosstalkFilenames = simSettings.crosstalkFilenames;
    zl = simSettings.zl;
    zs = simSettings.zs;
    
    % Display figure
    figure;
    hold on;

    % Plot main channel response
    PlotFreqResponse(channelFilename, zs, zl, 'Main Channel');

    % Plot crosstalk responses
    for i = 1:length(crosstalkFilenames)
        PlotFreqResponse(crosstalkFilenames{i}, zs, zl, ['Crosstalk ', num2str(i)]);
    end

    xlabel('Frequency (Hz)', 'Interpreter', 'latex');
    ylabel('Magnitude', 'Interpreter', 'latex');
    title('Frequency Responses of Channel and Crosstalk', 'Interpreter', 'latex');
    legend('Interpreter', 'latex');
    legend show;
    grid on;
    hold off;
end

function PlotFreqResponse(filename, zs, zl, label)
    sp = sparameters(filename);
    data = sp.Parameters;
    freq = sp.Frequencies;
    z0 = sp.Impedance;

    diffdata = s2sdd(data);
    diffz0 = 2 * z0;

    difftransfunc = s2tf(diffdata, diffz0, zs * 2, zl * 2);

    plot(freq*1.e-9,20*log10(abs(difftransfunc)),'DisplayName', label);
    xlim([0, max(freq)*1.e-9]);
end