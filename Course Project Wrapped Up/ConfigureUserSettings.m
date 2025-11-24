function [simSettings,simResources] = ConfigureUserSettings()
% Simulation Setting Configuration
    
    % Load S-parameter data (channel and crosstalk)
    % simSettings.channelFilename = 'DPO_10in_Meg7_THRU.s4p';
    % simSettings.crosstalkFilenames = {
    %                    'DPO_10in_Meg7_NENF1415.s4p',...
    %                    'DPO_10in_Meg7_NENF1112.s4p',...
    %                    'DPO_10in_Meg7_NENF89.s4p',...
    %                    'DPO_10in_Meg7_FENH1415.s4p',...
    %                    'DPO_10in_Meg7_FENH1112.s4p',...
    %                    'DPO_10in_Meg7_FENH89.s4p',...
    %                    'DPO_10in_Meg7_FENG1415.s4p',...
    %                    'DPO_10in_Meg7_FENG89.s4p'};

    % simSettings.channelFilename = 'DPO_12in_Meg7_THRU.s4p';
    % simSettings.crosstalkFilenames = {
    %                    'DPO_12in_Meg7_NENF1415.s4p',...
    %                    'DPO_12in_Meg7_NENF1112.s4p',...
    %                    'DPO_12in_Meg7_NENF89.s4p',...
    %                    'DPO_12in_Meg7_FENH1415.s4p',...
    %                    'DPO_12in_Meg7_FENH1112.s4p',...
    %                    'DPO_12in_Meg7_FENH89.s4p',...
    %                    'DPO_12in_Meg7_FENG1415.s4p',...
    %                    'DPO_12in_Meg7_FENG89.s4p'};

    % simSettings.channelFilename = 'DPO_14in_Meg7_THRU.s4p';
    % simSettings.crosstalkFilenames = {
    %                    'DPO_14in_Meg7_NENF1415.s4p',...
    %                    'DPO_14in_Meg7_NENF1112.s4p',...
    %                    'DPO_14in_Meg7_NENF89.s4p',...
    %                    'DPO_14in_Meg7_FENH1415.s4p',...
    %                    'DPO_14in_Meg7_FENH1112.s4p',...
    %                    'DPO_14in_Meg7_FENH89.s4p',...
    %                    'DPO_14in_Meg7_FENG1415.s4p',...
    %                    'DPO_14in_Meg7_FENG89.s4p'};

    simSettings.channelFilename = 'C2M__Z100_IL14_WC_BOR_H_L_H_THRU.s4p';
    simSettings.crosstalkFilenames = {
                       'C2M__Z100_IL14_WC_BOR_H_L_H_FEXT1.s4p',...
                       'C2M__Z100_IL14_WC_BOR_H_L_H_FEXT2.s4p',...
                       'C2M__Z100_IL14_WC_BOR_H_L_H_FEXT3.s4p',...
                       'C2M__Z100_IL14_WC_BOR_H_L_H_NEXT1.s4p',...
                       'C2M__Z100_IL14_WC_BOR_H_L_H_NEXT2.s4p',...
                       'C2M__Z100_IL14_WC_BOR_H_L_H_NEXT3.s4p',...
                       'C2M__Z100_IL14_WC_BOR_H_L_H_NEXT4.s4p'};

    % simSettings.channelFilename = 'TEC_Whisper27in_THRU_G14G15_07202016.s4p';
    % simSettings.crosstalkFilenames = {
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

    % Package Parameters
    simSettings.Cd = 100e-15;                                                   % Die-side capacitance (F) in π ladder model
    simSettings.Ls = 120e-12;                                                   % Series inductance (H) in π ladder model
    simSettings.Rs = 0;                                                         % Series resistance (Ohm) in π ladder model
    simSettings.Cb = 30e-15;                                                    % Bump-side capacitance (F) in π ladder model

    % Pre-cursors and post-cursors
    simSettings.trackPre  = 50;                                                 % Number of pre-cursors to track
    simSettings.trackPost = 100;                                                % Number of post-cursors to track
    simSettings.plotPre   = 5;                                                  % Number of pre-cursors to plot
    simSettings.plotPost  = 12;                                                 % Number of post-cursors to plot
    
    % Link settings
    simSettings.dataRate = 56 * 1e9;                                           % Data Rate (Hz)
    simSettings.pam = 4;                                                        % PAM levels (e.g. 2, 3, 4, 6, 8)
    simSettings.targetSER = 1e-4;                                               % Target SER
    simSettings.samplesPerSymb = 64;                                            % Samples per symbol
    simSettings.numSamples = 2^12;                                              % Number of samples to plot
    simSettings.numPlotPoints = 10000;                                          % Number of points to plot
    simSettings.jitterSampling = 0.01;                                          % Sampling jitter in UI (e.g., 0.01)
    simSettings.Tt = 0.4;                                                       % 20% - 80% rise/fall time in UI (e.g., 0.4)
    simSettings.txSNDR = 33;                                                    % In-band (Nyquist) TX SNDR in dB at 112 Gsym/s
    simSettings.rxPSD = 8e-9;                                                   % RX PSD (V^2/GHz)
    simSettings.diracSeparation = 0.02;                                         % Jitter pdf generation parameter in UI (e.g., 0.02)
    
    % Optimization Settings
    simSettings.opt.numTXFFE = 5;                                               % Number of precursor taps + 1 (main cursor)
    simSettings.opt.numRXFFE = 15;                                              % Number of RX FFE taps
    simSettings.opt.numRXDFE = 1;                                               % Number of RX DFE taps
    
    % Define parameters to sweep below in the format of (start : step : end)
    % or [val1, val2, val3, ...]

    % TX and RX termination resistance (Ohm)
    simSettings.sweep.txTerminationResistance = 50;
    simSettings.sweep.rxTerminationResistance = 50;
    
    % CTLE pole frequencies (Hz)
    % simSettings.sweep.ctle.fp1 = (20e9 : 500e7 : 30e9);
    % simSettings.sweep.ctle.fp2 = (500e6 : 100e6 : 900e6);
    % simSettings.sweep.ctle.fp3 = (50e9 : 1e9 : 56e9);
    simSettings.sweep.ctle.fp1 = (11.5e9 : 1000e7 : 21.5e9);
    simSettings.sweep.ctle.fp2 = (285e6 : 100e6 : 685e6);
    simSettings.sweep.ctle.fp3 = (28.6e9 : 2e9 : 36.6e9);
    
    % CTLE DC gains (dB)
    simSettings.sweep.ctle.g1 = (-20 : 5 : 0);
    simSettings.sweep.ctle.g2 = (-6 : 2 : 0);
    
    % RX front-end bandwidth limitation (Hz)
    % simSettings.sweep.rxBW = (25e9 : 1000e6 : 31e9);
    simSettings.sweep.rxBW = (25e9 : 1e9 : 31e9);

    % VGA gain (dB) and bandwidth (Hz)
    simSettings.sweep.vga.gain = (0 : 5 : 20);
    simSettings.sweep.vga.bw = 1e12;
    
    % Oversampling period
    simSettings.ts = 1/ simSettings.dataRate / simSettings.samplesPerSymb;

    % Simulation resources: initialize parallel pool
    if isempty(gcp('nocreate'))
        parpool;
    end
    
    % Simulation Resources: Do all of the rational fitting here
    [simSettings,simResources] = FitAllRationalModels(simSettings);
end


function [simSettings, simResources] = FitAllRationalModels(simSettings)
    cacheDir = 'rational_cache';
    if ~exist(cacheDir, 'dir')
        mkdir(cacheDir);  % Create a folder
    end

    %% === Main Channel ===
    baseName = matlab.lang.makeValidName(simSettings.channelFilename);
    cachePath = fullfile(cacheDir, [baseName '_main.mat']);

    if isfile(cachePath)
        fprintf('[INFO] Loading cached rational model for main channel: %s\n', simSettings.channelFilename);
        channel = load(cachePath, 'channel');
        simResources.channel = channel.channel;
        
        % Still need to load z0 separately
        S = sparameters(simSettings.channelFilename);
        simSettings.z0 = S.Impedance;
    else
        fprintf('[INFO] Fitting rational model for main channel: %s\n', simSettings.channelFilename);
        S = sparameters(simSettings.channelFilename);
        diffdata = s2sdd(S.Parameters);
        f = S.Frequencies;
        simSettings.z0 = S.Impedance;

        S_params = {
            squeeze(diffdata(1,1,:));
            squeeze(diffdata(1,2,:));
            squeeze(diffdata(2,1,:));
            squeeze(diffdata(2,2,:));
        };

        rS_params = cell(1,4);
        parfor k = 1:4
            rS_params{k} = rational(f, S_params{k}, 'Tolerance', -80);
        end

        channel.rS11 = rS_params{1};
        channel.rS12 = rS_params{2};
        channel.rS21 = rS_params{3};
        channel.rS22 = rS_params{4};

        simResources.channel = channel;
        save(cachePath, 'channel', '-v7.3');
        fprintf('[INFO] Saved main channel model: %s\n', cachePath);
    end

    %% === Crosstalk Channels ===
    numXT = numel(simSettings.crosstalkFilenames);
    filenames = simSettings.crosstalkFilenames;
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
            f = S.Frequencies;

            S_params = {
                squeeze(diffdata(1,1,:));
                squeeze(diffdata(1,2,:));
                squeeze(diffdata(2,1,:));
                squeeze(diffdata(2,2,:));
            };

            rS_params = cell(1,4);
            for k = 1:4
                rS_params{k} = rational(f, S_params{k}, 'Tolerance', -80);
            end

            localModel = struct();
            localModel.rS11 = rS_params{1};
            localModel.rS12 = rS_params{2};
            localModel.rS21 = rS_params{3};
            localModel.rS22 = rS_params{4};

            xtModels{i} = localModel;  % Safe to collect for post-parfor write
        end
    end

    % Save and assign to simResources
    simResources.crosstalk = repmat(struct('rS11', [], 'rS12', [], 'rS21', [], 'rS22', []), 1, numXT);
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
        simResources.crosstalk(i) = model;
    end
end