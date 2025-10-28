%-------------------------------------------------
% MAFTDSP Matlab Assignment 1 - Beyond the Basics
% 
% pitch shifting
% 
% Junting Wang 13/11/23
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
RA = 2;                                              % pitch shift ratio
RS = round(N/5);

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


% Calculate NFFT ----------------------------------------------------------
NFFT = 2^nextpow2(N);


% Initialise variables ----------------------------------------------------
phi_m = zeros(NFFT/2+1,1);
phi_m_plus_1 = zeros(NFFT/2+1,1);
theta_m = zeros(NFFT/2+1,1);
theta_m_plus_1 = zeros(NFFT/2+1,1);
omega_hat_k = (2 * pi * (0:NFFT/2)') / NFFT;


% Add analysis frame to the output vector ---------------------------------
for m = 1:NF-1
    xm_start = m * HA + 1;
    xm_end = xm_start + N-1;
    xm_win = win .* x(xm_start:xm_end);              % read in analysis frame and apply Hann window
    X = fft(xm_win, NFFT);                           % take the DFT of analysis frame
    X_half = X(1:NFFT/2+1);                          % keep only the first NFFT/2 + 1 bins
    Xmag = abs(X_half);                              % Separate into magnitude    
    Xang = angle(X_half);                            % Separate into phase
    phi_m_plus_1 = Xang;

    % Estimate the instantaneous frequencies
    delta_phi = ppa(phi_m_plus_1 - phi_m - omega_hat_k*HA); 
    phi_m = phi_m_plus_1;
    f_hat_m_plus_1 = delta_phi / HA;

    % Calculate the modified phases
    expected_phase = f_hat_m_plus_1 * HS;
    theta_m_plus_1 = ppa(theta_m + HS*omega_hat_k + expected_phase);
    theta_m = theta_m_plus_1;

    % Create new modified DFT frame Ym+1[k] of full-length NFFT
    Ym_plus_1 = Xmag .* exp(1j * theta_m_plus_1);
    Ym_plus_1_full = zeros(NFFT,1);
    Ym_plus_1_full(1:NFFT/2 + 1) = Ym_plus_1;
    Ym_plus_1_full(NFFT/2 + 2:end) = conj(Ym_plus_1(end-1:-1:2));  % Hermitian symmetry

    % Take the inverse DFT and resample the output frame
    ym = resample(real(ifft(Ym_plus_1_full, NFFT)),RA,1);

    % Truncate the frame from NFFT to N samples
    ym_truncated = ym(1:N);

    % Apply the synthesis window 
    ym_windowed = ym_truncated .* win;

    % Add the windowed synthesis frame into the output vector y
    y_start = m * RS + 1; 
    y_end = y_start + N -1;
    y(y_start:y_end) = y(y_start:y_end) + ym_windowed;
end


% listen to the output ----------------------------------------------------
soundsc(y,Fs);


% Create a ppa function ---------------------------------------------------
function wrapped_phase = ppa(phase)
wrapped_phase = mod(phase + pi, 2*pi) - pi;
end