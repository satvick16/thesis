function tapDFE = ExtractAlpha(simResult, simSettings)
% DFE tap extraction from cursor values

    cursorValues = simResult.cursorValues;
    numRXDFE = simSettings.opt.numRXDFE;

    % Extract max and its index
    [~, mainIdx] = max(cursorValues);

    % Extract DFE taps (CPU-side)
    tapDFE = cursorValues(mainIdx+1 : mainIdx+numRXDFE);
end
