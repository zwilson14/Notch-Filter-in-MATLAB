% Zach Wilson
% ENGE 330
% Final Project
% 12/3/16


% I didn't get time to create text prompts, 
% The first mouse click sets fc
% The second sets the high corner frequency, which sets bandwidth.

clear;

%%%%%%%%%%%%%%%%%%%%% Part A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[a,Fs] = audioread('testing.wav'); %Import the .wav file 
N = length(a);
tmax = (N-1)/Fs; 
t = linspace(0,tmax,N);
t = t'; %storing as column vector, to match x

% This FFT is for the Magnitude plot of a(t) -> A(f)
[AF, AF_freq] = fft330(a, Fs);

% This was for part A
noise = sin(600*2*pi.*t);
[NOISE, NOISE_f] = fft330((0.1*noise), Fs);
x = a + (0.1 .* noise);
[X, f] = fft330(x, Fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = ('Enter the number of harmonics:');
n_Max = inputdlg(prompt);
n_Max = str2num(n_Max{:})

%%% The Following is to set up the Fourier Series of 15.19 %%%

% These are the max values to set the parameters for the loops:
wo = 600 * 2 * pi; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change this used to be fc
c = 0; % This will hold the sawtooth function
A = 0.1; % Magnitude, Prevents sawtooth from over powering my voice

for n = 1:2:n_Max
    c = c + ((((2*A)/(pi*n))*sin(n*wo*t)) + (((-4*A)/(pi^2*n^2))*cos(n*wo*t)));
end

[C, CF] = fft330(c, Fs);

c_noise = c + a;
[C_Noise, C_f] = fft330(c_noise, Fs);


%%%%%%%%%%%%%%%%%%%%% Part C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

figure
plot(f, abs(C_Noise));
title('Frequency Domain of Inundated  nth Order Sawtooth With Voice Audio');
xlim([0, 1250]); % This defines the limits of the plot
% text_1 = text(f + 20, abs(X) 

% Set up mouse capture
[mouse_x, mouse_y] = ginput(2);
% These are used in Part B
fc = mouse_x(1);
fh = mouse_x(2); 
fl = (fc^2)/fh;
bw = abs(fh-fl);
%bw = abs((abs(fh) - abs(fc)) * 2);
% 1 means left button click, 2 is middle, 3 is right

% These set the values for the transfer function
w = 2*pi.*f;
wc = fc*2*pi;
s = j.*w;
z = bw/(2*fc); % zeta

order = n_Max; % This is set by the user
Hn = ones(size(a)); % creates a vector of ones

for n = 1:2:order
    fc_new = fc * n;
    wc_new = wc * n; %wc is the cutoff frequency determined by the user
    fl = ((fc_new)^2)/(fh*n);
    bw_new = abs(fh-fl);
    z_new = bw_new/(2*fc_new); % zeta
    
    H_new = (s.^2 + wc_new.^2) ./ (s.^2 + 2*z*wc_new.*s + wc_new.^2);
    Hn = Hn .* H_new;
end

Cleansed_Signal = (Hn .* (C_Noise));
[cleansed_signal, cleansed_freq] = ifft330(Cleansed_Signal, Fs); %time domain output

sound((cleansed_signal), Fs);

%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% These are the initial voice input in freq and time domain
figure
subplot(2, 1, 1);
plot(t, a);
title('Voice Input Without Noise - Time Domain');
xlabel('Time (s)');
ylabel('f(t)');

subplot(2, 1, 2);
plot(f, AF);
title('Voice Input Without Noise - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are the inundated signal in time and freq
figure
subplot(2, 1, 1);
plot(t, c_noise);
title('Time Domain of Inundated  nth Order Sawtooth With Voice Audio');
xlabel('Time (s)');
ylabel('f(t)');

subplot(2, 1, 2);
plot(f, C_Noise);
title('Frequency Domain of Inundated  nth Order Sawtooth With Voice Audio');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0, 1250]); % This defines the limits of the plot
ylim([-0.17, 0.125]); % This defines the limits of the plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are the Plot B spectral graphs, After the voice has been filtered
figure
subplot(2, 1, 1);
plot(t, cleansed_signal);
title('Filtered Signal - Time Domain');
xlabel('Time (s)');
ylabel('f(t)');

subplot(2, 1, 2);
plot(f, Cleansed_Signal);
title('Filtered Signal - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are the Plot B Sawtooth graphs
figure
subplot(2, 1, 1);
plot(t, (c)); % noise has its amplitude reduced
title('Sawtooth Noise - Time Domain');
xlabel('Time (s)');
ylabel('f(t)');

subplot(2, 1, 2);
plot(f, C);
title('Sawtooth Noise - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are the Plot B Transfer Function Plots
figure
plot(f, Hn);
title('nth Notch Filter - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




