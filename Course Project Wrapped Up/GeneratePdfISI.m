function simResult = GeneratePdfISI(simResult, simSettings, pltSettings)
% Generates the residual ISI PDF by convolving individual pre-cursor and post-cursor PDFs
% Inputs:
%   - simResult         : Simulation Results
%   - simSettings       : Simulation Settings
%   - pltSettings       : Plot Settings
%
% Output:
%   - pdfISI            : Generated ISI PDF (vector)

    cursorValues = simResult.cursorValues;
    pam = simSettings.pam;
    pdfRange = pltSettings.pdfRange;

    % Skip main cursor
    sigma = sum(cursorValues.^2) - max(cursorValues)^2;

    % Scale sigma by PAM levels and main cursor value
    sigma = sqrt(1/3 * (pam+1)/(pam-1) * sigma);
    pdfISI = normpdf(pdfRange, 0, sigma);

    simResult.pdfISI = pdfISI;
end

% function pdfISI = GeneratePdfISI(simResult, pltSettings)
% % Generates the residual ISI PDF by convolving individual pre-cursor and post-cursor PDFs using Dirac-style approximation
% %
% % INPUTS:
% %   simResult.cursorValues    : Vector of extracted cursor amplitudes
% %   pltSettings.pam           : Modulation order (e.g., 4, 8 for 4-PAM, 8-PAM)
% %   pltSettings.pdfRange      : X-axis for PDF evaluation (vector)
% %   pltSettings.pdfdx         : Bin width (scalar)
% %   pltSettings.trackPre      : Number of pre-cursor taps
% %   pltSettings.trackPost     : Number of post-cursor taps
% %
% % OUTPUT:
% %   pdfISI                    : Estimated ISI PDF (vector)
% 
%     cursorValues = simResult.cursorValues;
%     pam = pltSettings.pam;
%     pdfRange = pltSettings.pdfRange;
%     dx = pltSettings.pdfdx;
%     trackPre = pltSettings.trackPre;
%     trackPost = pltSettings.trackPost;
% 
%     % PAM symbol levels (normalized)
%     pamSymbols = linspace(-1, 1, pam);
% 
%     % Dirac-style single-cursor PDF builder
%     function pdf = cursor_pdf(h, pamSymbols, pdfRange, dx)
%         pdf = zeros(1, length(pdfRange));
%         for a = pamSymbols
%             x = h * a;
%             [~, idx] = min(abs(pdfRange - x));
%             pdf(idx) = pdf(idx) + 1 / (length(pamSymbols) * dx);
%         end
%     end
% 
%     % Construct ISI PDF by convolving each pre/post-cursor PDF
%     isi_indices = [1:trackPre, (trackPre+2):(trackPre+trackPost+1)];
%     pdfISI = cursor_pdf(cursorValues(isi_indices(1)), pamSymbols, pdfRange, dx);
% 
%     for k = 2:length(isi_indices)
%         h_k = cursorValues(isi_indices(k));
%         pdf_k = cursor_pdf(h_k, pamSymbols, pdfRange, dx);
%         pdfISI = conv(pdfISI, pdf_k, 'same') * dx;
%     end
% end
