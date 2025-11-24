function pltSettings = GeneratePltSettings(simResult)
% Plot Setting Configuration

    output = simResult.output;
    
    pdfRangeLimitScale = 10;
    pdfRangeResolution = 200001;

    pltSettings.pdfRange = linspace(-max(output)*pdfRangeLimitScale, max(output)*pdfRangeLimitScale, pdfRangeResolution);
    pltSettings.pdfdx = pltSettings.pdfRange(2) - pltSettings.pdfRange(1);
    
end