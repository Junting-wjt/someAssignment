%-------------------------------------------------
% MAFTDSP Matlab Assignment 2 - Part One
% 
% Comb filtering
% 
% Junting Wang 06/12/23
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
delay_time = 0.5;                                    % delay time [s]
g = 0.8;                                             % the 'strength' of the effect
M = round(Fs * delay_time);                          % integer number of samples in delay line


% Pre-allocate output vectors ---------------------------------------------
y_ff = zeros(size(x));                               % output of feedforward comb filter
y_fb = zeros(size(x));                               % output of feedback comb filter


% Pre-allocate delay-line buffer ------------------------------------------
dlinebuf = zeros(M,1);


% Perform feedforward comb filter in a for loop ---------------------------
for n = 1 : length(x)
    % apply the feedforword comb filter
    y_ff(n) = x(n) + g * dlinebuf(M);

    % update the delay line buffer
    dlinebuf(M) = x(n);
    dlinebuf = circshift(dlinebuf, 1);
end


% Perform feedback comb filter in a for loop ------------------------------
for n = 1 : length(x)
    % apply the feedback comb filter
    y_fb(n) = x(n) + g * dlinebuf(M);

    % update the delay line buffer
    dlinebuf(M) = y_fb(n);
    dlinebuf = circshift(dlinebuf, 1);
end


% Listen to the outputs ---------------------------------------------------
%soundsc(y_ff, Fs);
soundsc(y_fb, Fs);
