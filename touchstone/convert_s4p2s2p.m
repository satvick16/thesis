%% Convert 4-Port Single-Ended S-Parameters to 2-Port Differential S-Parameters
%
% -------------------------------------------------------------------------
% 1. PORT PAIRING FOR THIS FIXTURE
% -------------------------------------------------------------------------
% The original .s4p file contains the following header:
%
%   !! PORT DEFINITIONS:
%   !! Port 1 -----> TX Side G11     RX Side <----- Port 2
%   !! Port 3 -----> TX Side G12     RX Side <----- Port 4
%
% This means the true differential pairs are:
%
%   Differential Port 1 (TX):   Ports (1, 3)
%   Differential Port 2 (RX):   Ports (2, 4)
%
% MATLAB's s2sdd() **default option (option = 1)** uses the same pairing,
% so `s2sdd(s4p.Parameters)` automatically applies the correct mode
% conversion for this file.
%
% -------------------------------------------------------------------------
% 2. REFERENCE IMPEDANCE: 50 Ω (SE) → 100 Ω (DIFF)
% -------------------------------------------------------------------------
% The original S4P file is referenced to 50 Ω *per single-ended port*.
% A differential port consists of two 50-Ω single-ended conductors, so the
% matched differential impedance is:
%
%       Zdiff = Zse + Zse = 50 Ω + 50 Ω = 100 Ω
%
% When constructing the differential S2P file, we therefore set the
% reference impedance to 2*z0 = 100 Ω.

% Specify the S4P file name
input_file = 'DPO_10in_Meg7_THRU.s4p';

% Extract base name and create output .s2p filename
[folder, base, ~] = fileparts(input_file);
output_file = fullfile(folder, [base, '.s2p']);

% Load the 4-port single-ended S-parameter file
s4p = sparameters(input_file);

% Convert to 2-port
diffdata = s2sdd(s4p.Parameters);
freq = s4p.Frequencies;
z0 = s4p.Impedance;
s2p = sparameters(diffdata, freq, 2*z0);

% Write the output file
rfwrite(s2p, output_file);
fprintf('Converted "%s"  -->  "%s"\n', input_file, output_file);
