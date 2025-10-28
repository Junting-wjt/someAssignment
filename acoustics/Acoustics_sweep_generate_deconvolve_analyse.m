% TM UoE
% generate sweep, deconvolve, analyse

clear variables; clc ; %close all;
dbstop if error

%% GENERATE SWEEPS

sweep_length_s = 20;
freqRange = [20 20000];
fs = 48000;

%sweep_sig = generatesweep(freqRange(1), freqRange(2), sweep_length_s, fs, 0); % uncomment to generate sweep and inverse:

%% DECONVOLVE - Adjustable parameters

% initialInputGap = 0; % pre-delay before the detected onset for sweeps
% impulseLength_s = 10; % length of output impulse response in seconds
% impulseLength = impulseLength_s * fs; % length of impulse response in samples
% onsetThreshold = 0.5; % empirically chosen threshold used for truncation of IRs, chosen by visual observation of onset and noise floor amplitudes

clear irTrunc irAverage sweep sweep_aligned irRaw del del2 indexs indexMin;
%% Get on with the deconvolution now

% % read in generated sweeps:
% sweep_generated = audioread(['Sweep_20to20000_48000_',num2str(sweep_length_s),'s.wav']);
% inverseSweep = audioread(['InvSweep_20to20000_48000_',num2str(sweep_length_s),'s.wav']);

% read in generated sweeps:
sweep_generated = audioread('Sine_Wave_20s.wav');
inverseSweep = audioread('Sine_Wave_20s.wav');


% read in measured sweep(s):
%directory = strcat('C:\Users\86186\Desktop\Acoustics\acoustic_project_2023_12_05');
%[sweepMeas,fs] = audioread([directory,'Measured_pos1_sin20sec_with2people.wav']);
[sweepMeas,fs] = audioread("Sweep_IR_P2_20s.wav");
sweepMeas = sweepMeas(:,1); % just want first channel

% deconvolve:
irRaw = deconvolve(inverseSweep,sweepMeas);

%% Import previously measured IR

% [irRaw,fs] = audioread('/Users/mckenzt1/Documents/RT Dataset/rt dataset WAV/Wav Files/Meeting Room to Hallway/Source in Room/No Line of Sight/RIR_250cm.wav');
% [irRaw,fs] = audioread('/Users/mckenzt1/Downloads/st-georges-episcopal-church 2/stereo/st_georges_far.wav');
% [irRaw,fs] = audioread('IR01.wav');

%{
 [irRaw,fs] = audioread('IR_InFrontOfTable_Balloon_02.wav');
irRaw = irRaw(:,1); 
%}

impulseLength_s = length(irRaw)/fs;

%% ANALYSE

% spectrograms of generated sweep and inverse, measured sweep and deconvolved IR 
figure;
subplot(2,2,1); title('Original sweep');
spectrogram(sweep_generated,kaiser(256,5),220/2,512,fs,'yaxis');
ylim([0 20])

subplot(2,2,2); title('Inverse sweep');
spectrogram(inverseSweep,kaiser(256,5),220/2,512,fs,'yaxis');
ylim([0 20])

subplot(2,2,3); title('Measured sweep');
spectrogram(sweepMeas,kaiser(256,5),220/2,512,fs,'yaxis');
ylim([0 20])

subplot(2,2,4); title('Deconvolved Impulse Response');
spectrogram(irRaw',kaiser(256,5),220/2,512,fs,'yaxis');
ylim([0 20])

% plot frequency response:
figure; title('Frequency response of the Impulse Response');
f = fs/length(irRaw):fs/length(irRaw):fs; % Frequency vector for plotting
semilogx(f,20*log10(abs(fft(irRaw)))); % Plot the Fast Fourier transform
xlim([20 20000]); grid on; 
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

% spectrogram of the IR again but rotated and truncated
figure; title('Deconvolved Impulse Response');
% spectrogram(irRaw(1:impulseLength),kaiser(256,5),220/2,512,fs,'yaxis');
spectrogram(irRaw,kaiser(256,5),220/2,512,fs,'yaxis');
ylim([0 20])
view([58 28])

% plot time-domain amplitude response of the IR
t = (0:length(irRaw)-1)/fs';
figure;plot(t,irRaw); xlabel('Time (samples)'); ylabel('Amplitude');

%% Time-align and truncate
initialInputGap = 50; % pre-delay before the detected onset for sweeps
impulseLength_s = 10; % length of output impulse response in seconds
impulseLength = impulseLength_s * fs; % length of impulse response in samples
onsetThreshold = 0.1; % empirically chosen threshold used for truncation of IRs, chosen by visual observation of onset and noise floor amplitudes



irAbs = abs(irRaw); % absolute values of the IRs

% find when IR exceeds a set threshold amplitude
linearIndexes= find(irAbs>onsetThreshold);
indexs = linearIndexes(1);
indexMin = min(indexs);
% irTrunc = irRaw(indexMin - initialInputGap:(impulseLength+(indexMin - initialInputGap)-1));

irTrunc=irRaw; % uncomment to bypass truncation

% plot time-domain amplitude response of the truncated and aligned IR
figure;plot(irTrunc); xlabel('Time (samples)'); ylabel('Amplitude');

%% Plot Energy Decay Curve, linear regression and obtain RT60 value

% set parameters:
frequency = 2000;
rtFit = [-5 -25];

% bandpass filter
[B_decay, A_decay] = octdsgn(frequency, fs, 3);

% filter and schroeder integration to get EDC:
srir_input_filt = filter(B_decay,A_decay,irTrunc(:,1));
int_sch = 1/length(srir_input_filt) * cumtrapz(flip(srir_input_filt / max(abs(srir_input_filt))).^2);
edc = flip(int_sch(2:end));
edc = edc / max(edc); % normali  se
edc_dB = 10*log10(edc); % EDC in dB

% find values on EDC at which the decay has passed the values for the
% linear regression:
t1 = find(edc_dB <= rtFit(1),1,'first');
t2 = find(edc_dB <= rtFit(2),1,'first');

% get differences at these points, obtain regression:
x  = t2-t1; y = edc_dB(t2) - edc_dB(t1);
xy = y/x;
xvec = 1:1:length(edc_dB);
timevec = xvec./fs;
yvec = (1:1:length(edc_dB))*xy;

% align regression line to correct y axis value:
RT60_line = yvec+edc_dB(t1(1))-yvec(t1(1));
RT60_ind = find(RT60_line <= -60);
RT60 = timevec(RT60_ind(1)) % display rt60 value in commmand window

% plot EDC and regression:
figure;
plot(timevec,edc_dB,'k')
hold on;
ylim([-80 0]);
plot(timevec,RT60_line,'r--');
legend({'EDC','RT60'});
ylabel('Energy (dB)');xlabel('Time (s)');


%% Other possible metrics to calculate (you can do this yourselves!):
% Clarity, sound strength, direct-to-reverberant ratio (DRR), early decay
% time (EDT), multi-slope decays.
% [for spatial IRs] direction of arrival (DoA) of the direct sound
% and early reflections.  



