%-------------------------------------------------
% MAFTDSP Matlab Assignment 1 - Part One
% 
% A basic time stretcher
% 
% Junting Wang 11/11/23
%-------------------------------------------------


% Clear the command window, workspace and close all plots -----------------
clc;                                                 % clear the command window
clear;                                               % clear workspace
close all;                                           % close all plots


% Read in an input WAV file -----------------------------------------------
[stereoAudio, Fs] = audioread("mozart.wav");         % read in the audio file
x = (stereoAudio(:,1) + stereoAudio(:,2)) / 2;       % average left and right channels to mono


% Define analysis parameters ----------------------------------------------
window_time = 20;                                    % window time [ms]
O = 0.75;                                            % overlap factor O
Q = 1;                                               % time-stretch factor Q


% Determine variables -----------------------------------------------------                                                 
N = round((window_time / 1000) * Fs);                % window length N
HA = round(N - N * O);                               % analysis hop size HA
HS = round(Q * HA);                                  % synthesis stage hop size HS


% Generate a Hann window --------------------------------------------------
win = 0.5 * (1 - cos(2 * pi * (0:N - 1)'/ N));       % make the first value is zero 


% Determine the number of analysis frames ---------------------------------
for end_padding = 1:N-1
    x_padding = [zeros(N,1); x ;zeros(end_padding,1)];
    if mod(length(x_padding)-N,HA)==0
        end_padded = end_padding;                    % work out samples need to pad the end
    end
end
x=[zeros(N,1);x;zeros(end_padded,1)];
L = length(x);                                       % find the length of x
NF = 1+(L-N)/HA;


% Create an output vector y -----------------------------------------------
L_y = (NF-1)*HS + N;                                 % find the length of y
y = zeros(L_y,1);                                    % initially all zeros                               


% Add analysis frame into the output vector -------------------------------
for m = 1:NF-1
    % read in analysis frame
    xm_start = m * HA+1;                             
    xm_end = xm_start + N-1;

    % read in analysis frame and apply Hann window
    xm_winA = win.* x(xm_start:xm_end);             
    xm_winS = win.* xm_winA;

    % add to output vector
    y_start = m * HS+1;
    y_end = y_start + N - 1;
    y(y_start : y_end) = y(y_start : y_end) + xm_winS;
end


% listen to the output ----------------------------------------------------
soundsc(y,Fs);


% Verify exact reconstruction ---------------------------------------------
gain_factor = max(y)/max(x);
y = y ./ gain_factor;

if Q == 1
    t = (0:length(x)-1) * (1/Fs);
    figure(1);
    subplot(2,1,1);
    plot(t, x, 'r', t, y,'g');
    title('the Input and Output');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    subplot(2,1,2);
    error = y-x;
    plot(t, error,'b');
    title('the Error');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on; 
end


% plot spectrograms of the input and output vectors------------------------ 
figure(2);
subplot(2,1,1);
MA1_s2327978_Wang_myspec(x, Fs, N, O);
subplot(2,1,2);
MA1_s2327978_Wang_myspec(y, Fs, N, O);


%-------------------------------------------------
% For Q = 0.75, need a smaller N to maintain temporal resolution, keep
% the overlap value unchanged, and adjust HA through the values of N and O.
%
% For Q = 1.25, need a larger N to maintain frequency resolution, keep
% the overlap value unchanged, and adjust HA through the values of N and O.
%-------------------------------------------------