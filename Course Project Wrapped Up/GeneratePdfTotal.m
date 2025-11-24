function simResult = GeneratePdfTotal(simResult, simSettings, pltSettings, computeSER)
% Generates the total signal PDF by convolving all noise and jitter pdfs.
% Inputs:
%   - simResult             : Simulation Results
%   - simSettings           : Simulation Settings
%   - computeSER            : Flag that indicates SER estimation
% Output:
%   - pdfMainCursor         : Generated Main Cursor PDF (vector)
%   - pdfTotal              : Generated Total PDF (vector)
%   - SER                   : Estimated SER

    pam = simSettings.pam;
    pdfRange = pltSettings.pdfRange;
    pdfdx = pltSettings.pdfdx;
    cursorValues = simResult.cursorValues;
    pdfISI = simResult.pdfISI;
    pdfCrosstalk = simResult.pdfCrosstalk;
    pdfNoiseTX = simResult.pdfNoiseTX;
    pdfNoiseRX = simResult.pdfNoiseRX;
    pdfJitter = simResult.pdfJitter;
    

    % Extract the mainCursor
    mainCursor = max(cursorValues);

    % Generate the Dirac delta function at the cursor position
    pdfMainCursor = zeros(size(pdfRange));

    % Loop through each PAM level cursor
    for level = 0:pam-1
        % Calculate cursor level position (e.g., for PAM-2: -mainCursor and +mainCursor)
        inputLevel = -mainCursor + level * (2*mainCursor) / (pam-1);
        [~, idx] = min(abs(pdfRange - inputLevel));
        pdfMainCursor(idx) = 1 / pdfdx / pam;
    end
    
    % Convolve all noise PDFs together
    pdfShared = pdfISI;

    if ~isempty(pdfCrosstalk)
        pdfShared = conv(pdfShared, pdfCrosstalk, 'same') * pdfdx;
    end
    
    if ~isempty(pdfNoiseTX)
        pdfShared = conv(pdfShared, pdfNoiseTX, 'same') * pdfdx;
    end
    
    if ~isempty(pdfNoiseRX)
        pdfShared = conv(pdfShared, pdfNoiseRX, 'same') * pdfdx;
    end
    
    % === Combine per-level jitter PDFs with their Diracs ===
    pdfCond = cell(1, pam);
    pdfTotal = zeros(size(pdfRange));

    for level = 0:pam-1
        % Compute voltage position of this level
        inputLevel = -mainCursor + level * (2 * mainCursor) / (pam - 1);

        % Get corresponding jitter PDF
        jitterIdx = min(level+1, pam-level);    % symmetric mapping
        jitterPDF = pdfJitter{jitterIdx};

        % Convolve shared noise with level-specific jitter
        pdfLevel = conv(pdfShared, jitterPDF, 'same') * pdfdx;

        % Create Dirac at inputLevel
        pdfDirac = zeros(size(pdfRange));
        [~, idx] = min(abs(pdfRange - inputLevel));
        pdfDirac(idx) = 1 / pdfdx;

        % Convolve with Dirac
        pdfLevel = conv(pdfLevel, pdfDirac, 'same') * pdfdx;
        
        % Add to total pdf
        pdfTotal = pdfTotal + pdfLevel;
        
        % Store the conditional pdf
        pdfCond{level+1} = pdfLevel;
    end

    pdfTotal = pdfTotal / (sum(pdfTotal) * pdfdx);  % Normalize
    
    % For 2D pdf only: Estimate SER
    if computeSER
        eyeLevels = linspace(-mainCursor, mainCursor, pam);
        decisionBounds = [-Inf, 0.5 * (eyeLevels(1:end-1) + eyeLevels(2:end)), Inf];
        SER_per_symb = zeros(1,pam);
        
        for i = 1:pam
            cdf = cumsum(pdfCond{i}) * pdfdx;
            cdf = cdf / max(cdf);

            lowerB = decisionBounds(i);
            upperB = decisionBounds(i+1);

            idxL = find(pdfRange >= lowerB, 1, 'first');
            idxU = find(pdfRange >= upperB, 1, 'first');

            % Handle clipping
            if isempty(idxL), idxL = 1; end
            if isempty(idxU), idxU = length(cdf); end

            % Probability outside the decision region = left tail + right
            % tail
            SER_per_symb(i) = cdf(idxL) + (1 - cdf(idxU));
        end

        simResult.SER = mean(SER_per_symb);
    end

    simResult.pdfMainCursor = pdfMainCursor;
    simResult.pdfTotal = pdfTotal;
end
