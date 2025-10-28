%-------------------------------------------------
% MAFTDSP Matlab Assignment 2 - Part Two
% 
% Flanger Effect
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
delay_time = 0.001;                                  % delay time [s]
M0 = round(Fs * delay_time/2);                       % half the max number of samples of delay in the comb filter
f0 = 1;                                              % the frequency of LFO
g = 0.8;                                             % the 'strength' of the effect


% Pre-compute delay-profile vector M[n] -----------------------------------
M = round(M0 * (1+sin(2*pi*f0*(1:length(x))/Fs)));   % set a to one


% Pre-allocate the output vector ------------------------------------------
y = zeros(size(x));


% Pre-allocate a sufficiently large delay-line buffer ---------------------
dlinebuf = zeros(2 * M0 , 1);


% Perform the filtering operation in a for loop ---------------------------
for n = 1:length(x)
    % read the delay value and ensure samples not exceed the buffer limit
    delay_sample = dlinebuf(min(end,M(n)+1));

    % apply the effect
    y(n) = x(n) + g * delay_sample;

    % update the delay line buffer
    dlinebuf(end) = x(n);
    dlinebuf = circshift(dlinebuf, 1);
end


% Listen to the outputs ---------------------------------------------------
soundsc(y,Fs)
