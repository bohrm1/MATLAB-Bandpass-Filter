clear 
close all

%---------------------------filter paramters------------------------------%
Ft = 1000;                  % sampling freq

Fp_low = 410;               % lowpass passband freq 
Fs_low = 430;               % lowpass stopband freq 
delta_p_low = 0.01;         % lowpass passband ripple
delta_s_low = 0.1;          % lowpass stopband ripple

Fp_high = 60;               % highpass passband freq
Fs_high = 40;               % highpass stopband freq   
As_high = 50;               % highpass stopband attenuation
Ap_high = 2;                % highpass passband ripple


%------------------computing lowpass filter paramters---------------------%
% use firpmord to get order N_low and normalized edge frequencies Fo_low
[N_low, Fo_low, mo_low, W_low] = firpmord([Fp_low Fs_low],[1 0],[delta_p_low delta_s_low],Ft);
% use firmp function to generate LPF coeffs using result from firpmord
LPF_coeffs = firpm(N_low, Fo_low, mo_low, W_low);
% use dfilt.dffir to construct digital filter object from LPF_coeffs
LPF_filter = dfilt.dffir(LPF_coeffs);
% storing frequency response of LPF_filter in LPF_mag
LPF_mag = freqz(LPF_filter);


%------------------computing highpass filter paramters--------------------%
% computing ripples in linear scale 
delta_p_high = (10^(Ap_high/20)-1) / (10^(Ap_high/20)+1);
delta_s_high = (1+delta_p_high) / (10^(As_high/20));

% use firpmord to get order N_high and normalized edge frequencies Fo_high
[N_high, Fc_high, mo_high, W_high] = firpmord([Fs_high Fp_high],[0 1],[delta_s_high delta_p_high],Ft);
% use firmp function to generate HPF coeffs using result from firpmord
HPF_coeffs = firpm(N_high, Fc_high, mo_high, W_high);
% use dfilt.dffir to construct digital filter object from HPF_coeffs
HPF_filter = dfilt.dffir(HPF_coeffs);
% storing magnitude response of HPF_filter in HPF_mag
HPF_mag = freqz(HPF_filter);


%------------------computing bandpass filter responses--------------------%
% cascade LPF_filter and HPF_filter using dfilt.cascade() function
BPF_filter = dfilt.cascade(LPF_filter, HPF_filter);
% storing magnitude response of BPF_filter in BPF_mag
BPF_mag = freqz(BPF_filter);
% storing phase response of BPF_filter in BPF_phase
BPF_phase = phasez(BPF_filter);
% storing impulse response of BPF_filter in BPF_impulse
BPF_impulse = impz(BPF_filter);


%--------------------generating and filtering input-----------------------%
n = 0:Ft-1;                  % using 1000 samples for signal and FFT
k = 1:8192;        
k1 = k * (500/8192);           % normalizing n between 0 and 1 for FFT plots
x = cos(2*pi*(20/Ft)*n) + cos(2*pi*(40/Ft)*n) + cos(2*pi*(75/Ft)*n) + cos(2*pi*(100/Ft)*n) + cos(2*pi*(200/Ft)*n) + cos(2*pi*(300/Ft)*n) + cos(2*pi*(350/Ft)*n) + cos(2*pi*(440/Ft)*n) + cos(2*pi*(460/Ft)*n) + cos(2*pi*(475/Ft)*n);
x_fft = fft(x);              % FFT of input signal 

y = filter(BPF_filter, x);   % filtering x through BPF_filter and storing in y 
y_fft = fft(y);              % FFT of filtered output


%---------------------------generating plots------------------------------%
% figure 1: plot of input signal and its frequency content
figure(1); subplot(2,1,1); plot(n,x); title('x[n]'); xlim([0 500]), xlabel('Index'); ylabel('Amplitude');
           subplot(2,1,2); plot(n, abs(x_fft)); title('FFT x[n]'); xlim([0 500]), xlabel('Frequency (Hz)'); ylabel('Magnitude');

% figure2 : impulse response of BPF
figure(2); stem(BPF_impulse); title('BPF Impulse Response'), xlabel('Index'); ylabel('Amplitude');

% figure 3: BPF magnitude response (linear and dB) and BPF phase response
figure(3); subplot(3,1,1); plot(k1, abs(BPF_mag)); title('BPF Magnitude Response'); xlim([0 500]); xlabel('Frequency (Hz)'); ylabel('Magnitude');
           subplot(3,1,2); plot(k1, mag2db(abs(BPF_mag))); hold on; yline(-2, '--', '-2dB'); hold on; yline(-50, '--', '-50dB'); hold off; xlim([0 500]); ylim([-70 5]); title('BPF Magnitude Response (dB)'); xlabel('Frequency (Hz)'); ylabel('Magnitude');
           subplot(3,1,3); plot(k1, BPF_phase); title('BPF Phase Response'); xlim([0 500]); xlabel('Frequency (Hz)'); ylabel('Phase Delay');
           
% figure 4: plot of output signal and its frequency content           
figure(4); subplot(2,1,1); plot(n, y); title('y[n]'); xlim([0 500]), xlabel('Index'); ylabel('Amplitude');
           subplot(2,1,2); plot(n, abs(y_fft)); title('FFT y[n]'); xlim([0 500]), xlabel('Frequency (Hz)'); ylabel('Magnitude');
