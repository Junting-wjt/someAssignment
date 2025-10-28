%-------------------------------------------------
% MAFTDSP Matlab Assignment 1 - Part Three
% 
% My own spectrogram
% 
% Junting Wang 21/11/23
%-------------------------------------------------


function MA1_s2327978_Wang_myspec(x, Fs, N, overlap_factor)


% Determine variables -----------------------------------------------------
HA = round(N - N * overlap_factor);                         % analysis hop size HA
L = length(x);
NF = round(1+(L-N)/HA);                                     % the number of analysis frames


% Calculate NFFT ----------------------------------------------------------
NFFT = 2^nextpow2(N);


% Initialize the STFT matrix-----------------------------------------------
STFT_matrix = zeros(NFFT/2+1, NF);


% Generate a Hann window --------------------------------------------------
win = 0.5 * (1 - cos(2 * pi * (0:N - 1)'/ N));              % make the first value is zero 


% Store the DFT of frame into one column of STFT matrix-------------------- 
for m = 1:NF-1
    xm_start = (m-1) * HA+1;
    xm_end = xm_start + N-1;
    xm_win = win .* x(xm_start:xm_end);                     % read in analysis frame and apply Hann window
    X = fft(xm_win, NFFT);                                  % take the DFT of frame
    STFT_matrix(:,m) = X(1:NFFT/2+1);                       
end


% Plot the magnitude in decibels-------------------------------------------
STFT_dB = 20*log10(abs(STFT_matrix));
   

% Set frequency vector (zero to Nyquist frequency)-------------------------
f = (0:NFFT/2) * Fs / NFFT;


% Set time vector----------------------------------------------------------
t = (0:NF-1) * HA / Fs;
    

% Plot the spectrogram-----------------------------------------------------
imagesc(t, f, STFT_dB);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Calculate frame length and overlap factor in title
frame_length = N / Fs;
overlap = overlap_factor * 100;

% Set the title with frame length and overlap factor
str = sprintf('Frame Length: %.3f s, Overlap: %.1f%%',...
                  frame_length, overlap);
title(str);
    
colorbar;
colormap("hot");%set color mapping
    
% Find the maximum value in the STFT
maxdB = max(STFT_dB(:));

% Limit the color axis to a 60dB range from the highest decibels
caxis([maxdB-60 maxdB]);
end


