%% === Configure Simulation Settings ===
% Define channel and crosstalk aggressor file names
channelFilename = 'DPO_10in_Meg7_THRU.s4p';
crosstalkFilenames = {
                      'DPO_14in_Meg7_NENF1415.s4p',...
                      'DPO_14in_Meg7_NENF1112.s4p',...
                      'DPO_14in_Meg7_NENF89.s4p',...
                      'DPO_14in_Meg7_FENH1415.s4p',...
                      'DPO_14in_Meg7_FENH1112.s4p',...
                      'DPO_14in_Meg7_FENH89.s4p',...
                      'DPO_14in_Meg7_FENG1415.s4p',...
                      'DPO_14in_Meg7_FENG89.s4p'};

% Define number of pre-cursors and post-cursors to track and plot
trackPre = 50;                                                  % Number of pre-cursors to track
trackPost = 100;                                                % Number of post-cursors to track
plotPre = 5;                                                    % Number of pre-cursors to plot
plotPost = 12;                                                  % Number of post-cursors to plot

% Link settings
dataRate = 64 * 1e9;                                            % Data rate
pam = 4;                                                        % PAM levels (e.g. 2, 4, 6, 8)
SER = 1e-4;                                                     % Target SER
samplesPerSymb = 128;                                           % Samples per symbol
numSamples = 2^16;                                              % Number of samples to plot
numPlotPoints = 10000;                                          % Number of points to plot
jitterSampling = 0.01;                                          % Sampling jitter (UI)
txNoisePower = 1e-6;                                            % Transmitter Noise power
rxNoisePower = 1e-6;                                            % Receiver Noise power
txSNDR = 33;                                                    % TX SNDR (dB)
rxPSD = 8e-9;                                                   % RX PSD (V^2/GHz)
diracSeparation = 0.02;                                         % Jitter pdf generation parameter (e.g., 0.02)
addCrosstalk = true;                                            % Should add crosstalk
addNoise = true;                                                % Should add noise

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
fp1 = (11.5e9 : 1000e7 : 21.5e9);
fp2 = (285e6 : 100e6 : 685e6);
fp3 = (28.6e9 : 2e9 : 36.6e9);

% CTLE DC gains (dB)
g1 = (-20 : 5 : 0);
g2 = (-6 : 2 : 0);

% RX front-end bandwidth limitation (Hz)
% rxBW = (25e9 : 1000e6 : 31e9);
rxBW = (25e9 : 1e9 : 31e9);

% VGA gain (dB) and bandwidth (Hz)
vgaGain = (0 : 5 : 20);
vgaBW = 1e12;

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
fprintf('TX termination resistance sweep size: %d\n', length(txTerminationResistance));
fprintf('RX termination resistance sweep size: %d\n', length(rxTerminationResistance));
fprintf('CTLE f_p1 sweep size: %d\n', length(fp1));
fprintf('CTLE f_p2 sweep size: %d\n', length(fp2));
fprintf('CTLE f_p3 sweep size: %d\n', length(fp3));
fprintf('CTLE g1 sweep size: %d\n', length(g1));
fprintf('CTLE g2 sweep size: %d\n', length(g2));
fprintf('RX BW sweep size: %d\n', length(rxBW));
fprintf('VGA gain sweep size: %d\n', length(vgaGain));
fprintf('VGA bandwidth sweep size: %d\n', length(vgaBW));

% Flatten all sweep parameter combinations
[txZs, rxZs, g1s, g2s, fp1s, fp2s, fp3s, rxBWs, vgaGains, vgaBWs] = ndgrid( ...
    txTerminationResistance, rxTerminationResistance, ...
    g1, g2, fp1, fp2, fp3, ...
    rxBW, ...
    vgaGain, vgaBW);

totalIter = numel(txZs);
fprintf('Total iterations: %d\n', totalIter);

% Setup result tracking and progress reporting using DataQueue + handle class
dq = parallel.pool.DataQueue;
tracker = ProgressTracker(totalIter);
afterEach(dq, @(~) tracker.update());

%% === Simulation Begins ===
% === Generate Input (input is the time-domain input signal; t_in is the corresponding time vector)===
t_in = double((1:numSamples)') * ts;

% Generate symbol sequence: all 0s except one pulse at position 1
input_symbols = zeros(ceil(numSamples / samplesPerSymb), 1);
input_symbols(1) = 1;
input = kron(input_symbols, ones(samplesPerSymb, 1));

% Retrieve frequency response (main channel)
rS11 = channelModel.rS11;
rS12 = channelModel.rS12;
rS21 = channelModel.rS21;
rS22 = channelModel.rS22;

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
        % === Apply Channel ===
        fftFreqs = (0:numSamples-1)' / (ts * numSamples);
    
        % Evaluate frequency response
        S11 = squeeze(freqresp(rS11, fftFreqs));
        S12 = squeeze(freqresp(rS12, fftFreqs));
        S21 = squeeze(freqresp(rS21, fftFreqs));
        S22 = squeeze(freqresp(rS22, fftFreqs));
    
        % Compute H(f)
        Gamma1 = (2 * zs - 2 * z0) / (2 * zs + 2 * z0);
        Gamma2 = (2 * zl - 2 * z0) / (2 * zl + 2 * z0);
        DeltaS = S11 .* S22 - S12 .* S21;
        numerator = S21 .* (1 - Gamma1) .* (1 + Gamma2);
        denominator = 1 - S11 .* Gamma1 - S22 .* Gamma2 + Gamma1 * Gamma2 * DeltaS;
        Hf = numerator ./ denominator;
        
        % Apply channel to the input signal
        Xf = fft(input(:), numSamples);
        Yf = Hf(:) .* Xf;       % Channel output in frequency domain
    
        % === Apply CTLE === 
        if any([fp1, fp2, fp3] <= 0) || any(isnan([fp1, fp2, fp3])) ...
           || any(isnan([g1, g2])) || any(abs([g1, g2]) > 80)
            error('Iteration skipped due to invalid CTLE parameter values.');
        end
        
        % Compute CTLE transfer function H_ctle
        fz1 = fp1 * 10^(g1/20);
        fz2 = fp2 * 10^(g2/20);
        wz1 = 2 * pi * fz1;
        wz2 = 2 * pi * fz2;
    
        wp1 = 2 * pi * fp1;
        wp2 = 2 * pi * fp2;
        wp3 = 2 * pi * fp3;

        H_ctle = (1 + 1j*w/wz1) .* (1 + 1j*w/wz2) ./ ...
                 ((1 + 1j*w/wp1) .* (1 + 1j*w/wp2) .* (1 + 1j*w/wp3));
        
        % Apply CTLE in frequency domain
        Yf = Yf .* H_ctle;

        % === Apply Receiver Bandwidth Limitation ===
        H_bw = 1 ./ sqrt(1 + (fftFreqs(:) / rxBW).^(2 * 4));
        Yf = Yf .* H_bw;

        % === Apply VGA ===
        gain = 10^(vgaGain / 20);
        wc = 2 * pi * vgaBW;
        w = 2 * pi * fftFreqs(:);

        % Compute VGA transfer function: H(f) = gain / (j2πf + wc)
        H_vga = gain ./ sqrt(w.^2 + wc^2);

        % Apply VGA in frequency domain
        Yf = Yf .* H_vga;
        
        % Convert frequency domain signal to time domain
        output = real(ifft(Yf));

        % === Extract Cursors === 
        [~, mainIdx] = max(output);
        tapIndices = mainIdx + (-trackPre:trackPost) * samplesPerSymb;
        tapIndices = max(min(tapIndices, numel(output)), 1);
        cursorValues = output(tapIndices);      % [V]
        cursorTimes = t_ink(tapIndices) * 1e9;  % [ns]
        
        % === Adapt TX FFE ===
        
        % === Apply TX FFE ===
        
        % === Extract Cursors ===
    
        % === Extract Alpha ===
    
        % === Adapt RX FFE ===
    
        % === Apply RX FFE ===
    
        % === Extract Cursors and Define Decision Threshold ===
    
        % === Extract Alpha ===
    
        % === Apply DFE ===
    
        % === Extract Cursors ===
    
        % === Evaluate SNR ===
    
        % === Save Result ===

    catch err
        warning("Iteration %d failed: %s", idx, err.message);
        SNR_all(idx) = -Inf;
    end
end