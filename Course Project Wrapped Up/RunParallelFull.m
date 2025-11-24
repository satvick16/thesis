%% === Configure Simulation Settings === 
% Define channel and crosstalk aggressor file names
% channelFilename = 'DPO_10in_Meg7_THRU.s4p';
% crosstalkFilenames = {
%                       'DPO_14in_Meg7_NENF1415.s4p',...
%                       'DPO_14in_Meg7_NENF1112.s4p',...
%                       'DPO_14in_Meg7_NENF89.s4p',...
%                       'DPO_14in_Meg7_FENH1415.s4p',...
%                       'DPO_14in_Meg7_FENH1112.s4p',...
%                       'DPO_14in_Meg7_FENH89.s4p',...
%                       'DPO_14in_Meg7_FENG1415.s4p',...
%                       'DPO_14in_Meg7_FENG89.s4p'};

% channelFilename = 'DPO_12in_Meg7_THRU.s4p';
% crosstalkFilenames = {
%                    'DPO_12in_Meg7_NENF1415.s4p',...
%                    'DPO_12in_Meg7_NENF1112.s4p',...
%                    'DPO_12in_Meg7_NENF89.s4p',...
%                    'DPO_12in_Meg7_FENH1415.s4p',...
%                    'DPO_12in_Meg7_FENH1112.s4p',...
%                    'DPO_12in_Meg7_FENH89.s4p',...
%                    'DPO_12in_Meg7_FENG1415.s4p',...
%                    'DPO_12in_Meg7_FENG89.s4p'};

% channelFilename = 'DPO_14in_Meg7_THRU.s4p';
% crosstalkFilenames = {
%                    'DPO_14in_Meg7_NENF1415.s4p',...
%                    'DPO_14in_Meg7_NENF1112.s4p',...
%                    'DPO_14in_Meg7_NENF89.s4p',...
%                    'DPO_14in_Meg7_FENH1415.s4p',...
%                    'DPO_14in_Meg7_FENH1112.s4p',...
%                    'DPO_14in_Meg7_FENH89.s4p',...
%                    'DPO_14in_Meg7_FENG1415.s4p',...
%                    'DPO_14in_Meg7_FENG89.s4p'};

channelFilename = 'C2M__Z100_IL14_WC_BOR_H_L_H_THRU.s4p';
crosstalkFilenames = {
                   'C2M__Z100_IL14_WC_BOR_H_L_H_FEXT1.s4p',...
                   'C2M__Z100_IL14_WC_BOR_H_L_H_FEXT2.s4p',...
                   'C2M__Z100_IL14_WC_BOR_H_L_H_FEXT3.s4p',...
                   'C2M__Z100_IL14_WC_BOR_H_L_H_NEXT1.s4p',...
                   'C2M__Z100_IL14_WC_BOR_H_L_H_NEXT2.s4p',...
                   'C2M__Z100_IL14_WC_BOR_H_L_H_NEXT3.s4p',...
                   'C2M__Z100_IL14_WC_BOR_H_L_H_NEXT4.s4p'};

% channelFilename = 'TEC_Whisper27in_THRU_G14G15_07202016.s4p';
% crosstalkFilenames = {
%                    'TEC_Whisper27in_FEXT_F11F12_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_FEXT_F14F15_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_FEXT_F17F18_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_FEXT_G11G12_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_FEXT_G17G18_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_FEXT_H11H12_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_FEXT_H14H15_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_FEXT_H17H18_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_NEXT_F11F12_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_NEXT_F14F15_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_NEXT_F17F18_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_NEXT_G11G12_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_NEXT_G17G18_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_NEXT_H11H12_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_NEXT_H14H15_to_G14G15_07212016.s4p',...
%                    'TEC_Whisper27in_NEXT_H17H18_to_G14G15_07212016.s4p'};

% Package parameters
Cd = 100e-15;                                                   % Die-side capacitance (F) in π ladder model
Ls = 120e-12;                                                   % Series inductance (H) in π ladder model
Cb = 30e-15;                                                    % Bump-side capacitance (F) in π ladder model

% Define number of pre-cursors and post-cursors to track
trackPre = 50;                                                  % Number of pre-cursors to track
trackPost = 100;                                                % Number of post-cursors to track

% Link settings
dataRate = 56 * 1e9;                                            % Data rate
pam = 4;                                                        % PAM levels (e.g. 2, 4, 6, 8)
targetSER = 1e-4;                                               % Target SER
samplesPerSymb = 64;                                            % Samples per symbol
numSamples = 2^12;                                              % Number of samples to plot
jitterSampling = 0.01;                                          % Sampling jitter in UI (e.g., 0.01)
Tt = 0.4;                                                       % 20% - 80% rise/fall time in UI (e.g., 0.4)
txSNDR = 33;                                                    % In-band (Nyquist) TX SNDR in dB at 112 Gsym/s
rxPSD = 8e-9 * 1e-9;                                            % RX PSD (V^2/Hz)
diracSeparation = 0.02;                                         % Jitter pdf generation parameter in UI (e.g., 0.02)

% Optimization Settings
numTXFFE = 5;                                                   % Number of precursor taps + 1 (main cursor)
numRXFFE = 15;                                                  % Number of RX FFE taps
numRXDFE = 1;                                                   % Number of RX DFE taps

%% === Define parameters to sweep ===
% TX and RX termination resistance (Ohm)
txTerminationResistance = 50;
rxTerminationResistance = 50;

% CTLE pole frequencies (Hz)
% fp1 = (20e9 : 500e7 : 30e9);
% fp2 = (500e6 : 100e6 : 900e6);
% fp3 = (50e9 : 1e9 : 56e9);
fp1 = (20e9 : 2e9 : 30e9);
fp2 = (500e6 : 100e6 : 900e6);
fp3 = (50e9 : 2e9 : 56e9);

% CTLE DC gains (dB)
g1 = (-20 : 5 : 0);
g2 = (-6 : 2 : 0);

g1 = -10;
g2 = -4;
fp1 = dataRate / 2.5;
fp2 = dataRate / 80;
fp3 = dataRate;

% RX front-end bandwidth limitation (Hz)
% rxBW = (25e9 : 1000e6 : 31e9);
rxBW = (25e9 : 2e9 : 31e9);

rxBW = dataRate * 0.75;

% VGA gain (dB) and bandwidth (Hz)
vgaGain = (0 : 5 : 20);
vgaBW = 1e12;

vgaGain = 10;

%% === Simulation Setup (do NOT modify this section) ===
pulseWidth = 1 / dataRate;                                      % Input pulse width
ts = pulseWidth / samplesPerSymb;                               % Oversampling period

% Simulation resources: initialize parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

% Fit rational models for channels and crosstalk aggressors
cacheDir = 'rational_cache';
if ~exist(cacheDir, 'dir')
    mkdir(cacheDir);  % Create the folder if not exist
end

% Fit main channel
baseName = matlab.lang.makeValidName(channelFilename);
cachePath = fullfile(cacheDir, [baseName '_main.mat']);

if isfile(cachePath)
    fprintf('[INFO] Loading cached rational model for main channel: %s\n', channelFilename);
    channel = load(cachePath, 'channel');
    channelModel = channel.channel;
    
    % Still need to load z0 separately
    S = sparameters(channelFilename);
    z0 = S.Impedance;
else
    fprintf('[INFO] Fitting rational model for main channel: %s\n', channelFilename);
    S = sparameters(channelFilename);
    diffdata = s2sdd(S.Parameters);
    freq = S.Frequencies;
    z0 = S.Impedance;

    S_params = {
        squeeze(diffdata(1,1,:));
        squeeze(diffdata(1,2,:));
        squeeze(diffdata(2,1,:));
        squeeze(diffdata(2,2,:));
    };

    rS_params = cell(1,4);
    parfor k = 1:4
        rS_params{k} = rational(freq, S_params{k}, 'Tolerance', -80);
    end

    channel.rS11 = rS_params{1};
    channel.rS12 = rS_params{2};
    channel.rS21 = rS_params{3};
    channel.rS22 = rS_params{4};

    channelModel = channel;
    save(cachePath, 'channel', '-v7.3');
    fprintf('[INFO] Saved main channel model: %s\n', cachePath);
end

% === Fit Crosstalk Aggressors ===
numXT = numel(crosstalkFilenames);
filenames = crosstalkFilenames;
xtModels = cell(1, numXT);
xtCachePaths = cell(1, numXT);
xtNeedFit = false(1, numXT);

for i = 1:numXT
    baseName = matlab.lang.makeValidName(filenames{i});
    path = fullfile(cacheDir, [baseName '_xt.mat']);
    xtCachePaths{i} = path;
    if ~isfile(path)
        fprintf('[INFO] Crosstalk model NOT found. Will fit: %s\n', filenames{i});
        xtNeedFit(i) = true;
    else
        fprintf('[INFO] Crosstalk model cached. Will load: %s\n', filenames{i});
    end
end

% Fit only the models that are missing
parfor i = 1:numXT
    if xtNeedFit(i)
        file = filenames{i};
        S = sparameters(file);
        diffdata = s2sdd(S.Parameters);
        freq = S.Frequencies;

        S_params = {
            squeeze(diffdata(1,1,:));
            squeeze(diffdata(1,2,:));
            squeeze(diffdata(2,1,:));
            squeeze(diffdata(2,2,:));
        };

        rS_params = cell(1,4);
        for k = 1:4
            rS_params{k} = rational(freq, S_params{k}, 'Tolerance', -80);
        end

        localModel = struct();
        localModel.rS11 = rS_params{1};
        localModel.rS12 = rS_params{2};
        localModel.rS21 = rS_params{3};
        localModel.rS22 = rS_params{4};

        xtModels{i} = localModel;
    end
end

crosstalkModels = repmat(struct('rS11', [], 'rS12', [], 'rS21', [], 'rS22', []), 1, numXT);
for i = 1:numXT
    if xtNeedFit(i)
        model = xtModels{i};
        save(xtCachePaths{i}, 'model', '-v7.3');
        fprintf('[INFO] Saved crosstalk model: %s\n', filenames{i});
    else
        data = load(xtCachePaths{i}, 'model');
        model = data.model;
        fprintf('[INFO] Loaded crosstalk model: %s\n', filenames{i});
    end
    crosstalkModels(i) = model;
end

%% === Sanity Check: Display Sweep Sizes ===
fprintf('[INFO] TX termination resistance sweep size: %d\n', length(txTerminationResistance));
fprintf('[INFO] RX termination resistance sweep size: %d\n', length(rxTerminationResistance));
fprintf('[INFO] CTLE f_p1 sweep size: %d\n', length(fp1));
fprintf('[INFO] CTLE f_p2 sweep size: %d\n', length(fp2));
fprintf('[INFO] CTLE f_p3 sweep size: %d\n', length(fp3));
fprintf('[INFO] CTLE g1 sweep size: %d\n', length(g1));
fprintf('[INFO] CTLE g2 sweep size: %d\n', length(g2));
fprintf('[INFO] RX BW sweep size: %d\n', length(rxBW));
fprintf('[INFO] VGA gain sweep size: %d\n', length(vgaGain));
fprintf('[INFO] VGA bandwidth sweep size: %d\n', length(vgaBW));

% Flatten all sweep parameter combinations
[txZs, rxZs, g1s, g2s, fp1s, fp2s, fp3s, rxBWs, vgaGains, vgaBWs] = ndgrid( ...
    txTerminationResistance, rxTerminationResistance, ...
    g1, g2, fp1, fp2, fp3, ...
    rxBW, ...
    vgaGain, vgaBW);

totalIter = numel(txZs);
fprintf('\n[INFO] Total iterations: %d\n', totalIter);

% Setup result tracking and progress reporting using DataQueue + handle class
SNR_all = -Inf(totalIter, 1);
dq = parallel.pool.DataQueue;
tracker = ProgressTracker(totalIter);
afterEach(dq, @(~) tracker.update());

%% === Simulation Begins ===
% === Generate Input (input is the time-domain input signal; t_in is the corresponding time vector)===
outputTime = double((1:numSamples)') * ts;

% Generate symbol sequence: all 0s except one pulse at position 1
input_symbols = zeros(ceil(numSamples / samplesPerSymb), 1);
input_symbols(10) = 1;
input = kron(input_symbols, ones(samplesPerSymb, 1));
Xf = fft(input(:), numSamples);

% Retrieve frequency response (main channel)
rS11 = channelModel.rS11;
rS12 = channelModel.rS12;
rS21 = channelModel.rS21;
rS22 = channelModel.rS22;

% Generate FFT frequency vector (unshifted, matches fft() ordering)
Fs = 1 / ts;                  % sampling frequency
df = Fs / numSamples;         % frequency resolution
k  = (0:numSamples-1)';       % FFT bin indices
fftFreqs = k * df;
fftFreqs(k >= ceil(numSamples/2)) = fftFreqs(k >= ceil(numSamples/2)) - Fs; % Wrap bins above Nyquist to negative frequencies
w = 2 * pi * fftFreqs(:);

%% === Apply Driver Filter ===
a = 0.8 / (Tt / dataRate);
H_df = exp(-(pi * fftFreqs / a).^2);
Xf = H_df(:) .* Xf;

%% === Apply Package ===
H_tx_package = 1 ./ (1 + ((1j*w).^2) * (Ls * Cb));
H_rx_package = 1 ./ (1 + ((1j*w).^2) * (Ls * Cd));

Xf = H_tx_package(:) .* H_rx_package(:) .* Xf;

% Parallel Optimization
parfor idx = 1:totalIter
    % Unpack parameter combination
    zs = txZs(idx);
    zl = rxZs(idx);
    g1 = g1s(idx);
    g2 = g2s(idx);
    fp1 = fp1s(idx);
    fp2 = fp2s(idx);
    fp3 = fp3s(idx);
    rxBW = rxBWs(idx);
    vgaGain = vgaGains(idx);
    vgaBW = vgaBWs(idx);

    try
        %% === Apply Channel ===
        % Evaluate frequency response
        S11 = squeeze(freqresp(rS11, fftFreqs));
        S12 = squeeze(freqresp(rS12, fftFreqs));
        S21 = squeeze(freqresp(rS21, fftFreqs));
        S22 = squeeze(freqresp(rS22, fftFreqs));
    
        % Compute H_ch
        Gamma1 = (2 * zs - 2 * z0) / (2 * zs + 2 * z0);
        Gamma2 = (2 * zl - 2 * z0) / (2 * zl + 2 * z0);
        DeltaS = S11 .* S22 - S12 .* S21;
        numerator = S21 .* (1 - Gamma1) .* (1 + Gamma2);
        denominator = 1 - S11 .* Gamma1 - S22 .* Gamma2 + Gamma1 * Gamma2 * DeltaS;
        H_ch = numerator ./ denominator;
        
        % Apply channel to the input signal
        Yf = H_ch(:) .* Xf;       % Channel output in frequency domain

        %% === Apply Receiver Bandwidth Limitation ===
        x = fftFreqs(:) / rxBW;
        H_bw = 1 ./ (1 - 3.4142 * x.^2 + x.^4 + 1j * 2.6131 * (x - x.^3));
        Yf = Yf .* H_bw;

        %% === Apply CTLE === 
        if any([fp1, fp2, fp3] <= 0) || any(isnan([fp1, fp2, fp3])) ...
           || any(isnan([g1, g2])) || any(abs([g1, g2]) > 80)
            error('[ERROR] Iteration skipped due to invalid CTLE parameter values.');
        end
        
        % Compute CTLE transfer function H_ctle
        ctleGain = 10^((g1 + g2)/20);
        
        fz1 = fp1 * 10^(g1/20);
        fz2 = fp2 * 10^(g2/20);
        wz1 = 2 * pi * fz1;
        wz2 = 2 * pi * fz2;
        wp1 = 2 * pi * fp1;
        wp2 = 2 * pi * fp2;
        wp3 = 2 * pi * fp3;

        H_ctle = ctleGain .* (1 + 1j*w/wz1) .* (1 + 1j*w/wz2) ./ ...
                             ((1 + 1j*w/wp1) .* (1 + 1j*w/wp2) .* (1 + 1j*w/wp3));
        
        % Apply CTLE in frequency domain
        Yf = Yf .* H_ctle;

        %% === Apply VGA ===
        wc = 2 * pi * vgaBW;

        % Compute VGA transfer function: H_vga(f) = gain / (j2πf + wc)
        H_vga = 10^(vgaGain / 20) * (wc ./ (1j * w + wc));

        % Apply VGA in frequency domain
        Yf = Yf .* H_vga;
        
        % Convert frequency domain signal to time domain
        output = real(ifft(Yf));

        %% === Extract Cursors === 
        [~, mainIdx] = max(output);
        tapIndices = mainIdx + (-trackPre:trackPost) * samplesPerSymb;
        validMask = tapIndices >= 1 & tapIndices <= numel(output);
        tapIndices = tapIndices(validMask);
        cursorValues = output(tapIndices);              % [V]

        %% === Adapt TX FFE ===
        % Extract segment with main + (N-1) pre-cursors
        [~, mainIdx] = max(cursorValues);
        firstIdx = max(1, mainIdx - numTXFFE + 1);
        h = cursorValues(firstIdx:mainIdx);     % h = [h_{-N+1}, ..., h0]
        h = h(:).';                             % Ensure row vector
    
        % Construct pulse response matrix H (N × (2N-1))
        N = numTXFFE;
        H = zeros(N, 2*N - 1);
        for i = 1:N
            H(i, i:i + N - 1) = h;
        end
    
        % Construct target pulse response g*
        g_star = zeros(1, size(H, 2));
        g_star(end) = h(end);  % Desired main cursor value at the end
    
        % Compute optimal FFE taps
        tx_ffe_taps = g_star * H' / (H * H');
    
        % Normalize to maintain consistent TX swing
        tx_ffe_taps = tx_ffe_taps / sum(abs(tx_ffe_taps));
        
        %% === Apply TX FFE ===
        % Upsample FFE taps to match oversampling rate
        tapsFFE_upsampled = zeros(1, (length(tx_ffe_taps) - 1) * samplesPerSymb + 1);
        tapsFFE_upsampled(1:samplesPerSymb:end) = tx_ffe_taps;
        
        % Apply FFE to main signal
        output = conv(output, tapsFFE_upsampled, 'same');

        %% === Extract Cursors ===
        [~, mainIdx] = max(output);
        tapIndices = mainIdx + (-trackPre:trackPost) * samplesPerSymb;
        validMask = tapIndices >= 1 & tapIndices <= numel(output);
        tapIndices = tapIndices(validMask);
        cursorValues = output(tapIndices);              % [V]

        %% === Extract Alpha ===
        [~, mainIdx] = max(cursorValues);
        tapsDFE = cursorValues(mainIdx+1 : mainIdx+numRXDFE);

        %% === Adapt RX FFE ===
        % Compute RX noise power
        H_total = abs(H_bw .* H_ctle .* H_vga).^2;      % Total frequency Response
        df = 1 / (ts * numSamples);                     % Frequency resolution
        sigma_n2 = rxPSD * sum(H_total) * df;           % Total shaped noise power
        
        % Initialization
        minMMSE = Inf;
        optSamplingPoint = 0;
        w_opt = zeros(numRXFFE, 1);
        
        % Construct matrix X (Toeplitz-like)
        X = zeros(numRXFFE, numRXFFE + length(cursorValues) - 1);
        for row = 1:numRXFFE
            X(row, row : row + length(cursorValues) - 1) = cursorValues;
        end
        y_length = size(X, 2);
        
        % Sweep over candidate sampling points
        for i = 1:y_length
            % Build ideal target response
            y = zeros(y_length, 1);
            y(i) = 1;
    
            % Add post-cursors for DFE co-optimization
            for j = 1:length(tapsDFE)
                index = i + j;
                if index <= y_length
                    y(index) = tapsDFE(j);
                end
            end
    
            % Build system A w = b
            A = X * X.';
            b = X * y;
    
            % Check conditioning and apply regularization if needed
            if rcond(A) < 1e-6
                reg_strength = sigma_n2;
                A = A + reg_strength * eye(size(A));
            end
    
            % Solve system
            try
                w_i = A \ b;
            catch
                continue;  % Skip this sampling point if solving fails
            end
    
            % Compute MMSE
            y_hat = X.' * w_i;
            mmse = norm(y_hat - y)^2 + sigma_n2 * norm(w_i)^2;
    
            % Track best result
            if mmse < minMMSE
                minMMSE = mmse;
                w_opt = w_i;
                optSamplingPoint = i;
            end
        end
    
        % Normalize
        if norm(w_opt) > 0
            rx_ffe_taps = w_opt / norm(w_opt);
        else
            rx_ffe_taps = w_opt;
        end

        %% === Apply RX FFE ===
        % Upsample FFE taps to match oversampling rate
        tapsFFE_upsampled = zeros(1, (length(rx_ffe_taps) - 1) * samplesPerSymb + 1);
        tapsFFE_upsampled(1:samplesPerSymb:end) = rx_ffe_taps;
        
        % Apply FFE to main signal
        output = conv(output, tapsFFE_upsampled, 'same');

        %% === Extract Cursors and Define Decision Threshold ===
        [mainValue, mainIdx] = max(output);
        tapIndices = mainIdx + (-trackPre:trackPost) * samplesPerSymb;
        validMask = tapIndices >= 1 & tapIndices <= numel(output);
        tapIndices = tapIndices(validMask);
        cursorValues = output(tapIndices);              % [V]
        % cursorTimes = outputTime(tapIndices) * 1e9;     % [ns]

        % Set DFE decision threshold
        thresh = mainValue / 2;

        %% === Extract Alpha ===
        [~, mainIdx] = max(cursorValues);
        tapsDFE = cursorValues(mainIdx+1 : mainIdx+numRXDFE);

        %% === Apply DFE ===
        [~, mainIdx] = max(output);
        outputDFE = output;
        detected_symbols = zeros(size(output));
        n_start = mod(mainIdx - 1, samplesPerSymb) + 1;
        
        for n = n_start:samplesPerSymb:length(output)
            % Decision device with threshold
            if outputDFE(n) > thresh
                detected_symbols(n) = 1;
            else
                detected_symbols(n) = -1;
            end
    
            % DFE in discrete time
            for k = 1:length(tapsDFE)
                delayed_idx = n + (k-1/2) * samplesPerSymb;  % Delay by k*UI
                if delayed_idx <= length(output) && detected_symbols(n) == 1
                    % Spread correction effect over symbol interval
                    idx_range = delayed_idx : min(delayed_idx + samplesPerSymb - 1, length(output));
                    outputDFE(idx_range) = outputDFE(idx_range) - tapsDFE(k) * detected_symbols(n);
                end
            end
        end
        output = outputDFE;

        %% === Extract Cursors ===
        [mainValue, mainIdx] = max(output);
        tapIndices = mainIdx + (-trackPre:trackPost) * samplesPerSymb;
        validMask = tapIndices >= 1 & tapIndices <= numel(output);
        tapIndices = tapIndices(validMask);
        cursorValues = output(tapIndices);              % [V]

        %% === Evaluate SNR ===
        % Compute Noise Variance
        P_signal = (1/3) * (pam + 1) / (pam - 1);
        P_noise_tx = P_signal / (10^(txSNDR / 10));     % TX noise power at transmitter
        txPSD = P_noise_tx / 56e9;
        
        w_ffe = 2 * pi * fftFreqs(:) / dataRate;
        E_tx = exp(-1j * w_ffe * (0:(numel(tx_ffe_taps)-1)));
        E_rx = exp(-1j * w_ffe * (0:(numel(rx_ffe_taps)-1)));

        H_tx_ffe = E_tx * tx_ffe_taps(:);
        H_rx_ffe = E_rx * rx_ffe_taps;

        H_tx = H_ch(:) .* H_bw(:) .* H_ctle(:) .* H_vga(:);     % Composite TX noise transfer function
        G_tx = sum(abs(H_tx_ffe).^2 .* abs(H_tx).^2 .* abs(H_rx_ffe).^2) * df;

        H_rx = H_bw(:) .* H_ctle(:) .* H_vga(:);                % Composite RX noise transfer function
        G_rx = sum(abs(H_rx).^2 .* abs(H_rx_ffe).^2) * df;

        TXNoiseVar = txPSD * G_tx;
        RXNoiseVar = rxPSD * G_rx;

        pamVarianceNormFactor = 1/3 * (pam + 1) / (pam - 1);
        signal_power = mainValue^2 * pamVarianceNormFactor;
        residual_power = sum(cursorValues.^2) * pamVarianceNormFactor - signal_power;
        SNR = 10 * log10(signal_power / (residual_power + TXNoiseVar + RXNoiseVar));

        %% === Save Result ===
        SNR_all(idx) = SNR;

    catch err
        warning("[WARNING] Iteration %d failed: %s", idx, err.message);
        SNR_all(idx) = -Inf;
    end

    % Report progress
    send(dq, []);
end

%% === Report Best Result ===
[maxSNR, bestIdx] = max(SNR_all);

fprintf('\n === Max SNR achieved: %.4f dB ===\n', maxSNR);
fprintf('\tzs = %d Ohm\n', txZs(bestIdx));
fprintf('\tzl = %d Ohm\n', rxZs(bestIdx));
fprintf('\tg1 = %d dB\n', g1s(bestIdx));
fprintf('\tg2 = %d dB\n', g2s(bestIdx));
fprintf('\tfp1 = %d Hz\n', fp1s(bestIdx));
fprintf('\tfp2 = %d Hz\n', fp2s(bestIdx));
fprintf('\tfp3 = %d Hz\n', fp3s(bestIdx));
fprintf('\tReceiver Bandwidth = %d Hz\n', rxBWs(bestIdx));
fprintf('\tVGA Gain = %d dB\n', vgaGains(bestIdx));
fprintf('\tVGA Bandwidth = %d Hz\n', vgaBWs(bestIdx));