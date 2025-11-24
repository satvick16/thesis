function [cursorTimes, cursorValues] = ExtractCursors(simResult, simSettings)
% Extracts cursor tap times and values from the output response.
%
% INPUTS:
%   simResult.output           - Output waveform
%   simResult.outputTime       - Time vector (in seconds)
%   simSettings.samplesPerSymb - Oversampling ratio
%   simSettings.pulseWidth     - Symbol period (in seconds)
%   simSettings.trackPre       - Number of precursor taps
%   simSettings.trackPost      - Number of postcursor taps
%
% OUTPUTS:
%   cursorTimes   - Sampling times (in nanoseconds)
%   cursorValues  - Corresponding amplitudes

    output = simResult.output;
    outputTime = simResult.outputTime;  % [s]
    sps = simSettings.samplesPerSymb;
    trackPre = simSettings.trackPre;
    trackPost = simSettings.trackPost;

    % Find main cursor (maximum)
    [~, mainIdx] = max(output);

    % Cursor index offsets
    tapOffsets = (-trackPre:trackPost) * sps;
    tapIndices = mainIdx + tapOffsets;

    % Keep only valid indices
    validMask = tapIndices >= 1 & tapIndices <= numel(output);
    tapIndices = tapIndices(validMask);

    % Extract values
    cursorValues = output(tapIndices);
    cursorTimes = outputTime(tapIndices) * 1e9;  % [ns]
end

% function [cursorTimes, cursorValues] = ExtractCursors(simResult, simSettings)
% % Extracts cursor tap times and values from the output response.
% %
% % INPUTS:
% %   simResult.output           - Output waveform
% %   simResult.outputTime       - Time vector (in seconds)
% %   simSettings.samplesPerSymb - Oversampling ratio
% %   simSettings.pulseWidth     - Symbol period (in seconds)
% %   simSettings.trackPre       - Number of precursor taps
% %   simSettings.trackPost      - Number of postcursor taps
% %
% % OUTPUTS:
% %   cursorTimes   - Sampling times (in nanoseconds)
% %   cursorValues  - Corresponding amplitudes
% 
%     output = simResult.output;
%     outputTime = simResult.outputTime;  % Assumed to be in seconds
%     sps = simSettings.samplesPerSymb;
% 
%     % Find main cursor (maximum)
%     [mainVal, mainIdx] = max(output);
%     mainTime = outputTime(mainIdx) * 1e9;  % Convert to ns
% 
%     % Initialize
%     cursorTimes = nan(1, simSettings.trackPre + 1 + simSettings.trackPost);
%     cursorValues = nan(1, simSettings.trackPre + 1 + simSettings.trackPost);
% 
%     % Main cursor
%     center = simSettings.trackPre + 1;
%     cursorTimes(center) = mainTime;
%     cursorValues(center) = mainVal;
% 
%     % Precursors
%     for k = 1:simSettings.trackPre
%         idx = mainIdx - k * sps;
%         if idx < 1, continue; end
%         cursorTimes(center - k) = outputTime(idx) * 1e9;
%         cursorValues(center - k) = output(idx);
%     end
% 
%     % Postcursors
%     for k = 1:simSettings.trackPost
%         idx = mainIdx + k * sps;
%         if idx > length(output), continue; end
%         cursorTimes(center + k) = outputTime(idx) * 1e9;
%         cursorValues(center + k) = output(idx);
%     end
% end

% function [cursorTimes, cursorValues] = ExtractCursors(simResult, simSettings)
% % GPU-compatible version of ExtractCursors
% 
%     output = simResult.output;
%     outputTime = simResult.outputTime;
%     sps = simSettings.samplesPerSymb;
% 
%     % Step 1: find main peak (needs CPU for indexing)
%     [mainVal, mainIdx] = max(gather(output));  % extract index from GPU
%     mainTime = gather(outputTime(mainIdx)) * 1e9;
% 
%     % Step 2: initialize cursor arrays (CPU-side)
%     numCursors = simSettings.trackPre + 1 + simSettings.trackPost;
%     cursorTimes = nan(1, numCursors);
%     cursorValues = nan(1, numCursors);
% 
%     center = simSettings.trackPre + 1;
%     cursorTimes(center) = mainTime;
%     cursorValues(center) = gather(mainVal);
% 
%     % Step 3: precursors
%     for k = 1:simSettings.trackPre
%         idx = mainIdx - k * sps;
%         if idx < 1, continue; end
%         cursorTimes(center - k) = gather(outputTime(idx)) * 1e9;
%         cursorValues(center - k) = gather(output(idx));
%     end
% 
%     % Step 4: postcursors
%     for k = 1:simSettings.trackPost
%         idx = mainIdx + k * sps;
%         if idx > length(output), continue; end
%         cursorTimes(center + k) = gather(outputTime(idx)) * 1e9;
%         cursorValues(center + k) = gather(output(idx));
%     end
% end