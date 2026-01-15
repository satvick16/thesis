% 1. Define the filename
filename = 'channel_two_port.s2p'; 

% 2. Import the S-parameter data into an S-parameter object
% Make sure 'your_file_name.s2p' is in your current MATLAB directory or specify the full path.
S = sparameters(filename);

% 3. Plot all the S-parameters
% This will generate a standard plot (e.g., magnitude in dB vs. frequency).
rfplot(S);

% Optional: Add a grid to the plot
grid on;