% Define the file name
filename = 'eric.wav';

% Read the audio file
[audioSignal, Fs] = audioread(filename);

% Define the length of the audio signal
L = length(audioSignal);

% Define the time vector for plotting
t = (0:L-1) / Fs;

% Compute the FFT of the audio signal
Y = fft(audioSignal);

% Compute the frequency vector
f = (0:L-1) * (Fs / L);
f_shifted = (-L/2:L/2-1) * (Fs / L); % Frequency vector for plotting

% Shift zero frequency component to the center
Y_shifted = fftshift(Y);

% Plot the spectrum of the original audio signal
figure;
plot(f_shifted, abs(Y_shifted));
title('Spectrum of the Original Audio Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-Fs/2 Fs/2]); % Set x-axis limits to show the full frequency range
grid on;

% Define the cutoff frequency for the low-pass filter
cutoff_frequency = 4e3; % 4 kHz

% Create an ideal low-pass filter
filter = zeros(size(Y));
% Set filter to 1 for frequencies <= cutoff_frequency
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

% Plot the filtered signal in the time and frequency domains
% Time-Domain Plot
figure;
subplot(2,1,1);
plot(t, audioSignal_filtered);
title('Filtered Audio Signal in Time Domain (BW = 4 kHz)');
xlabel('Time (s)');
ylabel('Amplitude');

% Frequency-Domain Plot (Zoomed to ±5 kHz)
subplot(2,1,2);
plot(f_shifted, abs(Y_filtered_shifted));
title('Filtered Audio Signal in Frequency Domain (BW = 4 kHz)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-5e3 5e3]); % Focus on the range of interest around the cutoff frequency
grid on;

% Ensure the filtered signal is in the correct range for audio playback
audioSignal_filtered_normalized = audioSignal_filtered / max(abs(audioSignal_filtered));

% Play the filtered audio signal
disp('Playing the filtered signal...');
sound(audioSignal_filtered_normalized, Fs);

% Parameters for modulation
Fc = 100e3;            % Carrier frequency (100 kHz)
Fs_new = 5 * Fc;       % New sampling frequency (5 times the carrier frequency)

% resampling method
resample_factor = Fs_new / Fs;
t_resampled = (0:round(length(audioSignal_filtered_normalized) * resample_factor)-1) / Fs_new;
audioSignal_resampled = interp1((0:length(audioSignal_filtered_normalized)-1) / Fs, audioSignal_filtered_normalized, t_resampled, 'linear', 0);

% Generate the time vector for the resampled signal
t_resampled = (0:length(audioSignal_resampled)-1) / Fs_new;
% DSB-TC Modulation (Standard AM Modulation)
modulation_index = 0.5; % Given modulation index
carrier_signal = cos(2 * pi * Fc * t_resampled); % Carrier signal
message_signal = audioSignal_resampled; % Message signal

% Add the carrier component to the message signal based on modulation index
DSB_TC_signal = (1 + modulation_index * message_signal) .* carrier_signal; % DSB-TC modulation (standard AM)

% DSB-SC Modulation (Suppressed Carrier)
DSB_SC_signal = message_signal .* carrier_signal; % DSB-SC modulation

% Plot the modulated signals in the time domain
figure;
subplot(2,1,1);
plot(t_resampled, DSB_TC_signal);
title('DSB-TC Modulated Signal (Time Domain)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_resampled, DSB_SC_signal);
title('DSB-SC Modulated Signal (Time Domain)');
xlabel('Time (s)');
ylabel('Amplitude');

% Compute and plot the spectrum of the modulated signals in the frequency domain
N = length(t_resampled); % Number of samples for FFT
DSB_TC_spectrum = fftshift(fft(DSB_TC_signal, N)); % FFT and shift for DSB-TC
DSB_SC_spectrum = fftshift(fft(DSB_SC_signal, N)); % FFT and shift for DSB-SC
f_resampled = (-N/2:N/2-1) * (Fs_new / N); % Frequency vector for the new sampling rate

% Plot the spectra
figure;
subplot(2,1,1);
plot(f_resampled, abs(DSB_TC_spectrum));
title('DSB-TC Modulated Signal (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(2,1,2);
plot(f_resampled, abs(DSB_SC_spectrum));
title('DSB-SC Modulated Signal (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

xlim([-2*Fc 2*Fc]); % Zoom in around the carrier frequencies
grid on;

% Envelope Detection for DSB-TC
envelope_DSB_TC = abs(hilbert(DSB_TC_signal));

% Envelope Detection for DSB-SC
envelope_DSB_SC = abs(hilbert(DSB_SC_signal));

% Plot the envelopes
figure;
subplot(2,1,1);
plot(t_resampled, envelope_DSB_TC);
title('Envelope of DSB-TC Modulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_resampled, envelope_DSB_SC);
title('Envelope of DSB-SC Modulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Downsample the envelope of DSB-TC signal
received_DSB_TC = resample(envelope_DSB_TC, Fs, Fs_new);

% Downsample the envelope of DSB-SC signal
received_DSB_SC = resample(envelope_DSB_SC, Fs, Fs_new);

% Normalize the received DSB-TC signal
received_DSB_TC_normalized = received_DSB_TC / max(abs(received_DSB_TC));

% Normalize the received DSB-SC signal
received_DSB_SC_normalized = received_DSB_SC / max(abs(received_DSB_SC));

% Time vector for plotting received signals
t_received = (0:length(received_DSB_TC)-1) / Fs;

figure;
subplot(2,1,1);
plot(t_received, received_DSB_TC_normalized);
title('Received DSB-TC Signal in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_received, received_DSB_SC_normalized);
title('Received DSB-SC Signal in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

% Pause to allow the filtered signal to finish playing
pause(length(received_DSB_TC_normalized) / Fs + 2);  % Adding 2 seconds pause

% Play the received DSB-TC signal
disp('Playing the received DSB-TC signal...');
sound(received_DSB_TC_normalized, Fs);

% Pause to allow the DSB-TC signal to finish playing
pause(length(received_DSB_TC_normalized) / Fs + 2);  % Adding 2 seconds pause

% Play the received DSB-SC signal
disp('Playing the received DSB-SC signal...');
sound(received_DSB_SC_normalized, Fs);

% After playing both signals, you can observe the following:
% The DSB-TC signal should be accurately demodulated by the envelope detector, as the carrier helps retrieve the original signal.
% The DSB-SC signal, however, may not be demodulated correctly, since the envelope detector does not perform well without a carrier signal. The output might sound distorted or incorrect.

% This observation shows that the envelope detector is suitable for DSB-TC modulation but not for DSB-SC modulation.


% Define SNR values
snr_values = [0, 10, 30];

% Initialize storage for received signals with different SNR
received_signals_snr = cell(1, length(snr_values));

% Coherent Detection for DSB-SC with different SNR values
for i = 1:length(snr_values)
    % Add noise to the DSB-SC signal with specified SNR additive wide gauss
    % noise
    noisy_DSB_SC_signal = awgn(DSB_SC_signal, snr_values(i), 'measured');

    % Coherent detection (demodulation)
    coherent_demodulated_signal = noisy_DSB_SC_signal .* cos(2 * pi * Fc * t_resampled);

    % Low-pass filter the demodulated signal to remove high-frequency components
    [b, a] = but
    ter(5, cutoff_frequency/(Fs_new/2));  % 5th-order Butterworth filter
    received_signal_snr = filtfilt(b, a, coherent_demodulated_signal);
    
    % Downsample the received signal back to the original sampling rate
    received_signal_snr = resample(received_signal_snr, Fs, Fs_new);
    
    % Normalize the received signal for playback
    received_signal_snr_normalized = received_signal_snr / max(abs(received_signal_snr));
    
    % Store the received signal
    received_signals_snr{i} = received_signal_snr_normalized;
    
    % Plot the received signal in the time domain
    figure;
    subplot(2,1,1);
    plot(t_received, received_signal_snr_normalized);
    title(['Received Coherent Signal with SNR = ', num2str(snr_values(i)), ' dB (Time Domain)']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % Plot the spectrum of the received signal in the frequency domain
    received_spectrum_snr = fftshift(fft(received_signal_snr, N));
    subplot(2,1,2);
    plot(f_resampled, abs(received_spectrum_snr));
    title(['Received Coherent Signal Spectrum with SNR = ', num2str(snr_values(i)), ' dB (Frequency Domain)']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([-2*Fc 2*Fc]);  % Zoom in around the carrier frequencies
    grid on;
    
    pause(length(received_signal_snr_normalized) / Fs + 2);  % Wait for the sound to finish playing
    % Play the received signal
    disp(['Playing the received DSB-SC signal with SNR = ', num2str(snr_values(i)), ' dB...','with FC = 100K Hz']);
    sound(received_signal_snr_normalized, Fs);
end

pause(length(received_signal_snr_normalized) / Fs + 2);  % Wait for the sound to finish playing

% Define the frequency error
Fc_error = 100.1e3;  % New carrier frequency with error of 100.1 kHz

% Define the SNR values for the test
snr_values_error = [0, 10, 30];  % SNR values for 0 dB, 10 dB, and 30 dB

for i = 1:length(snr_values_error)
    % Add noise to the DSB-SC signal with specified SNR
    noisy_DSB_SC_signal_error = awgn(DSB_SC_signal, snr_values_error(i), 'measured');

    % Coherent detection with frequency error (demodulation)
    coherent_demodulated_signal_error = noisy_DSB_SC_signal_error .* cos(2 * pi * Fc_error * t_resampled);

    % Low-pass filter the demodulated signal to remove high-frequency components
    received_signal_error_snr = filtfilt(b, a, coherent_demodulated_signal_error);
    
    % Downsample the received signal back to the original sampling rate
    received_signal_error_snr = resample(received_signal_error_snr, Fs, Fs_new);
    
    % Normalize the received signal for playback
    received_signal_error_snr_normalized = received_signal_error_snr / max(abs(received_signal_error_snr));
    
    % Store the received signal with frequency error
    received_signals_error_snr{i} = received_signal_error_snr_normalized;
    
    % Plot the received signal with frequency error in the time domain
    figure;
    subplot(2,1,1);
    plot(t_received, received_signal_error_snr_normalized);
    title(['Received DSB-SC Signal with Frequency Error = 100.1 kHz and SNR = ', num2str(snr_values_error(i)), ' dB (Time Domain)']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % Plot the spectrum of the received signal with frequency error in the frequency domain
    received_spectrum_error_snr = fftshift(fft(received_signal_error_snr, N));
    subplot(2,1,2);
    plot(f_resampled, abs(received_spectrum_error_snr));
    title(['Received DSB-SC Signal with Frequency Error = 100.1 kHz and SNR = ', num2str(snr_values_error(i)), ' dB (Frequency Domain)']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([-2*Fc 2*Fc]);  % Zoom in around the carrier frequencies
    grid on;
    
    pause(length(received_signal_error_snr_normalized) / Fs + 2);  % Wait for the sound to finish playing
    % Play the received signal with frequency error
    disp(['Playing the received DSB-SC signal with Frequency Error = 100.1 kHz and SNR = ', num2str(snr_values_error(i)), ' dB...']);
    sound(received_signal_error_snr_normalized, Fs);
end

pause(length(received_signal_error_snr_normalized) / Fs + 2);  % Wait for the sound to finish playing

% Calculate the frequency error
frequency_error = abs(Fc - Fc_error);

% Display the frequency error
disp(['The frequency error is: ', num2str(frequency_error), ' Hz']);
fprintf('The phenomenon caused by the frequency error is known as Carrier Frequency Offset (CFO).');

% Define the phase error
phase_error_degrees = 20; % Phase error in degrees
phase_error_radians = deg2rad(phase_error_degrees); % Convert to radians

% Define the SNR values for the test
snr_values_phase_error = [0, 10, 30];  % SNR values for 0 dB, 10 dB, and 30 dB

for i = 1:length(snr_values_phase_error)
    % Add noise to the DSB-SC signal with specified SNR
    noisy_DSB_SC_signal_phase_error = awgn(DSB_SC_signal, snr_values_phase_error(i), 'measured');

    % Coherent detection with phase error (demodulation)
    % Introduce phase error into the cosine term
    coherent_demodulated_signal_phase_error = noisy_DSB_SC_signal_phase_error .* cos(2 * pi * Fc * t_resampled + phase_error_radians);

    % Low-pass filter the demodulated signal to remove high-frequency components
    received_signal_phase_error_snr = filtfilt(b, a, coherent_demodulated_signal_phase_error);
    
    % Downsample the received signal back to the original sampling rate
    received_signal_phase_error_snr = resample(received_signal_phase_error_snr, Fs, Fs_new);
    
    % Normalize the received signal for playback
    received_signal_phase_error_snr_normalized = received_signal_phase_error_snr / max(abs(received_signal_phase_error_snr));
    
    % Store the received signal with phase error
    received_signals_phase_error_snr{i} = received_signal_phase_error_snr_normalized;
    
    % Plot the received signal with phase error in the time domain
    figure;
    subplot(2,1,1);
    plot(t_received, received_signal_phase_error_snr_normalized);
    title(['Received DSB-SC Signal with Phase Error = 20 degrees and SNR = ', num2str(snr_values_phase_error(i)), ' dB (Time Domain)']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % Plot the spectrum of the received signal with phase error in the frequency domain
    received_spectrum_phase_error_snr = fftshift(fft(received_signal_phase_error_snr, N));
    subplot(2,1,2);
    plot(f_resampled, abs(received_spectrum_phase_error_snr));
    title(['Received DSB-SC Signal with Phase Error = 20 degrees and SNR = ', num2str(snr_values_phase_error(i)), ' dB (Frequency Domain)']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([-2*Fc 2*Fc]);  % Zoom in around the carrier frequencies
    grid on;
    
    pause(length(received_signal_phase_error_snr_normalized) / Fs + 2);  % Wait for the sound to finish playing
    % Play the received signal with phase error
    disp(['Playing the received DSB-SC signal with Phase Error = 20 degrees and SNR = ', num2str(snr_values_phase_error(i)), ' dB...']);
    sound(received_signal_phase_error_snr_normalized, Fs);
end

pause(length(received_signal_phase_error_snr_normalized) / Fs + 2);  % Wait for the sound to finish playing

% Display the phase error information
disp(['The phase error is: ', num2str(phase_error_degrees), ' degrees']);
fprintf('The phenomenon caused by the phase error can affect the signal’s phase alignment and decoding accuracy.');
