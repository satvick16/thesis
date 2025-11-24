% Load the 4-port single-ended S-parameter file
s4p = sparameters('C2M__Z100_IL14_WC_BOR_H_L_H_THRU.s4p');

% Extract S-parameter matrix and frequencies
S = s4p.Parameters;     % 4x4xN
freq = s4p.Frequencies; % Nx1
Z0 = s4p.Impedance;

% Define the single-ended to differential transformation matrix
% Differential mode voltage = (V+ - V-) / sqrt(2)
% Common mode voltage = (V+ + V-) / sqrt(2)
T = 1/sqrt(2) * [ 1 -1  0  0;
                  0  0  1 -1;
                  1  1  0  0;
                  0  0  1  1 ];

% Apply the mixed-mode transformation for each frequency point
Nf = length(freq);
Smm = zeros(4,4,Nf);
for k = 1:Nf
    Smm(:,:,k) = T * S(:,:,k) * inv(T);
end

% Extract only the differential-to-differential (Sdd) 2x2 submatrix
Sdd = Smm(1:2,1:2,:);

% Create new 2-port S-parameter object
s2p = sparameters(Sdd, freq, Z0);

% Save as Touchstone .s2p file
rfwrite(s2p, 'C2M__Z100_IL14_WC_BOR_H_L_H_THRU_single_ended.s2p');
