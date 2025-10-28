%-------------------------------------------------
% MAFTDSP Matlab Assignment 2 - Part Three
% 
% Chorus Effect
% 
% Junting Wang 08/12/23
%-------------------------------------------------


% Clear the command window, workspace and close all plots -----------------
clc;                                                 % clear the command window
clear;                                               % clear workspace
close all;                                           % close all plots


% Read in an input WAV file -----------------------------------------------
[x, Fs] = audioread("Cath_cut.wav");                 % read in the audio file
[r_number, c_number] = size(x);
if c_number == 2
    x = (x(:,1) + x(:,2)) / 2;                       % average left and right channels to mono
end


% Define parameters -------------------------------------------------------
delay_time1 = 0.025;                                 % delay line time for LFO1 [s]
swing_time1 = 0.008;                                 % swing range time for LFO1 [s]
f1 = 1;                                              % LFO frequency for LFO1
g1 = 0.9;                                            % effect strength for LFO1
delay_time2 = 0.035;                                 % delay line time for LFO2 [s]
swing_time2 = 0.010;                                 % swing range time for LFO2 [s]
f2 = 2;                                              % LFO frequency for LFO2
g2 = 0.8;                                            % effect strength for LFO2


% Convert effect depths and swing ranges to samples -----------------------
P01 = round(delay_time1 * Fs);                       % effect depth for LFO1 in samples
D1 = round(swing_time1 * Fs);                        % effect oscillation for LFO1 in samples
P02 = round(delay_time2 * Fs);                       % effect depth for LFO2 in samples
D2 = round(swing_time2 * Fs);                        % effect oscillation for LFO2 in samples


% Pre-compute delay-profile vectors M1 and M2 -----------------------------
M1 = round(P01 + D1 * sin(2 * pi * f1 * (1:length(x)) / Fs));
M2 = round(P02 + D2 * sin(2 * pi * f2 * (1:length(x)) / Fs));


% Pre-allocate the output vector ------------------------------------------
y = zeros(size(x));


% Pre-allocate a sufficiently large delay-line buffer ---------------------
max_delay = max([M1, M2]);
dlinebuf = zeros(max_delay,1);


% Perform the chorus filtering operation in a for loop --------------------
for n = 1:length(x)
    % read from vectors M1 and M2
    current_delay1 = M1(n);
    current_delay2 = M2(n);

    % ensure index not exceed the buffer limit and no invalid index values
    index1 = max(1, min(length(dlinebuf), current_delay1));
    index2 = max(1, min(length(dlinebuf), current_delay2));

    % apply the effect
    y(n) = x(n) + g1 * dlinebuf(index1) + g2 * dlinebuf(index2);

    % update the delay line buffer
    dlinebuf(end) = x(n);
    dlinebuf = circshift(dlinebuf, 1);
end


% Listen to the outputs ---------------------------------------------------
sound(y, Fs);
