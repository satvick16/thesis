% Manual conversion from 4-port differential S-parameters to 2-port single-ended S-parameters
% No RF Toolbox needed

inputFile  = 'C2M__Z100_IL14_WC_BOR_H_L_H_THRU.s4p';   % Input differential S-parameter file
outputFile = 'channel_two_port.s2p';   % Output single-ended file

% --- Read .s4p file ---
fid = fopen(inputFile);
lines = {};
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if startsWith(line, '!') || isempty(line)
        continue;
    end
    lines{end+1} = line; %#ok<SAGROW>
end
fclose(fid);

% Extract frequency and S-parameters
header = split(lines{1});
freq_unit = header{2}; % e.g. 'GHz'

data = [];
for i = 2:length(lines)
    nums = str2double(split(lines{i}));
    data = [data; nums]; %#ok<AGROW>
end

freq = data(:,1);
S = zeros(4,4,length(freq));  % 4x4 S-matrix

% Parse each line of magnitude/angle data
% Touchstone order: S11 S12 S13 S14 S21 S22 ... S44
for k = 1:length(freq)
    idx = 2;
    for r = 1:4
        for c = 1:4
            mag = data(k,idx); ang = data(k,idx+1);
            S(r,c,k) = mag .* exp(1j*ang*pi/180);
            idx = idx + 2;
        end
    end
end

% --- Perform mixed-mode to single-ended conversion ---
M = (1/sqrt(2)) * [1  0  1  0;
                   -1 0  1  0;
                    0  1  0  1;
                    0 -1  0  1];

S2 = zeros(2,2,length(freq));

for k = 1:length(freq)
    S2(:,:,k) = M * S(:,:,k) * inv(M);
    % Take the 2x2 differential block
    S2(:,:,k) = S2(1:2,1:2);
end

% --- Write output .s2p file ---
fid = fopen(outputFile, 'w');
fprintf(fid, '# %s S MA R 50\n', freq_unit);
for k = 1:length(freq)
    fprintf(fid, '%e ', freq(k));
    for r = 1:2
        for c = 1:2
            mag = abs(S2(r,c,k));
            ang = angle(S2(r,c,k))*180/pi;
            fprintf(fid, '%e %e ', mag, ang);
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);

fprintf('Converted %s → %s\n', inputFile, outputFile);
