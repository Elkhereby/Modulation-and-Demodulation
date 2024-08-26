% Define the file name
filename = 'eric.wav';

% Read the audio file
[audioSignal, Fs_original] = audioread(filename);

% Define the length of the audio signal
L = length(audioSignal);

% Define the time vector for plotting
t = (0:L-1) / Fs_original;

% Compute the FFT of the audio signal
Y = fft(audioSignal);

% Compute the frequency vector
f = (0:L-1) * (Fs_original / L);
f_shifted = (-L/2:L/2-1) * (Fs_original / L); % Frequency vector for plotting

% Shift zero frequency component to the center
Y_shifted = fftshift(Y);

% Plot the spectrum of the original audio signal
figure;
plot(f_shifted, abs(Y_shifted));
title('Spectrum of the Original Audio Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-Fs_original/2 Fs_original/2]); % Set x-axis limits to show the full frequency range
grid on;

% Define the cutoff frequency for the low-pass filter
cutoff_frequency = 4e3; % 4 kHz

% Create an ideal low-pass filter
filter = zeros(size(Y));
filter(abs(f_shifted) <= cutoff_frequency) = 1; 

% Apply the filter in the frequency domain
Y_filtered = Y_shifted .* filter;

% Perform inverse FFT to get the filtered time-domain signal
audioSignal_filtered = real(ifft(ifftshift(Y_filtered)));

% Plot original and filtered signals
figure;
subplot(2,1,1);
plot(t, audioSignal);
title('Original Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, audioSignal_filtered);
title('Filtered Audio Signal (Low-pass filter at 4 kHz)');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the spectrum of the filtered signal (zoomed to show ±5 kHz)
Y_filtered_final = fft(audioSignal_filtered);
Y_filtered_shifted = fftshift(Y_filtered_final);
figure;
plot(f_shifted, abs(Y_filtered_shifted));
title('Spectrum of the Filtered Audio Signal (BW = 4 kHz)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-5e3 5e3]); % Zoom in to the ±5 kHz range to observe the filter effect
grid on;

% Ensure the filtered signal is in the correct range for audio playback
audioSignal_filtered_normalized = audioSignal_filtered / max(abs(audioSignal_filtered));

% Play the filtered audio signal
disp('Playing the filtered signal...');
sound(audioSignal_filtered_normalized, Fs_original);

% Pause to allow the filtered signal to finish playing
pause(length(audioSignal_filtered_normalized) / Fs_original + 2); % Adding 2 seconds extra buffer time

%% FM Modulation Section
% Clear any previous variables that might interfere
clear y F t1;

% Calculate modulation index kf
kf = 75;  % Adjust as necessary

% Define carrier frequency and amplitude
fc = 100000;  % Carrier frequency 100 kHz
Ac = 1;       % Carrier amplitude

% Resample the filtered signal for modulation
y_resampled = resample(audioSignal_filtered, 5*fc, Fs_original);
Fs_resampled = 5*fc;  % Update sampling frequency after resampling

% Calculate the time vector for the resampled signal
t1 = linspace(0, length(y_resampled) / Fs_resampled, length(y_resampled))';

% Perform FM modulation
%cumsum is summtion equivalent to integeration
FM_signal = Ac * cos(2*pi*fc*t1 + 2*pi*kf*cumsum(y_resampled)/Fs_resampled);

% Fourier transform of the FM signal
L_FM = length(FM_signal);
F_FM = fftshift(fft(FM_signal));
f_FM = Fs_resampled/2 * linspace(-1, 1, L_FM);

% Plot the FM modulated signal in the frequency domain
figure;
plot(f_FM, abs(F_FM));
title('FM Modulation in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% FM Demodulation Section
% Clear previous variables that may interfere
clear envelopeFM;

% Apply a differentiator (discriminator)
d_y = diff(FM_signal);
d_y = [0; d_y];  % Pad with zero to maintain the length

% Envelope detection using Hilbert transform
envelopeFM = abs(hilbert(d_y)) - mean(abs(hilbert(d_y)));

% Plot the demodulated signal in the time domain
figure;
plot(t1, envelopeFM);
title('Demodulated FM in Time Domain using Envelope Detector');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-2*10^-4 2*10^-4]);

% Resample the demodulated signal back to the original sampling frequency
envelopeFM_resampled = resample(envelopeFM, Fs_original, Fs_resampled);

% Normalize the signal for playback
envelopeFM_normalized = envelopeFM_resampled / max(abs(envelopeFM_resampled));

% Play the demodulated signal
disp('Playing the demodulated signal...');
sound(500 * abs(envelopeFM_normalized), Fs_original);  % Adjust the volume factor as needed
