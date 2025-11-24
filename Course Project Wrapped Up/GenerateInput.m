function simResult = GenerateInput(simSettings, should_plot)
% GENERATE_INPUT Generates a single-bit pulse input signal for time-domain simulation
%
% INPUTS:
%   simSettings.dataRate        - Data rate (Hz)
%   simSettings.samplesPerSymb  - Number of samples per symbol
%   simSettings.numSamples      - Total number of samples
%   should_plot                 - If true, plot the generated signal
%
% OUTPUTS:
%   output                      - Time-domain signal (vector)
%   Yf                          - Frequency-domain signal (vector)

    dataRate = simSettings.dataRate;
    samplesPerSymb = simSettings.samplesPerSymb;
    numSamples = simSettings.numSamples;
    ts = simSettings.ts;

    % Generate FFT frequency vector (unshifted, matches fft() ordering)
    Fs = 1 / ts;                  % sampling frequency
    df = Fs / numSamples;         % frequency resolution
    k  = (0:numSamples-1)';       % FFT bin indices
    fftFreqs = k * df;

    % Wrap bins above Nyquist to negative frequencies
    fftFreqs(k >= ceil(numSamples/2)) = fftFreqs(k >= ceil(numSamples/2)) - Fs;
    simResult.fftFreqs = fftFreqs;
    
    % Time vector
    outputTime = double((1:numSamples)') * ts;

    % Generate symbol sequence: all 0s except one pulse at position 1
    input_symbols = zeros(ceil(numSamples / samplesPerSymb), 1);
    input_symbols(10) = 1;

    % Expand to continuous-time input signal
    input = kron(input_symbols, ones(samplesPerSymb, 1));
    
    simResult.Yf = fft(input, numSamples);
    simResult.input = input;    % This will be used as the time-domain excitation for crosstalk simulation
    simResult.output = input;
    simResult.outputTime = outputTime;

    % Optional plotting
    if should_plot
        pulseWidth = 1 / dataRate;
        figure;
        plot(outputTime(1:samplesPerSymb * 5) * 1e9, input(1:samplesPerSymb * 5), 'LineWidth', 2);
        title([num2str(dataRate * 1e-9), ' Gbps Input Signal'], 'FontSize', 12);
        ylabel('Input signal');
        xlabel('Time (ns)');
        axis([-inf, inf, -0.5, 1.5]);
        hold on;

        % Symbol cursor lines
        symbol_times = ((0:4) + 0.5) * pulseWidth * 1e9;
        for t = symbol_times
            xline(t, 'Color', [0.5, 0.5, 0.5]);
        end
        hold off;
    end
end