[audio_waveform, fs_audio] = audioread ('eric.wav'); % read the audio file

% Fourier transform
audio_spectrum = fftshift (fft(audio_waveform) );
f = fs_audio/2 * linspace(-1, 1, length (audio_waveform) ) ;

% Plot the original signal in the frequency domain
figure;
subplot (2,2,1);
plot(f, abs (audio_spectrum) / length (audio_waveform) ) ;
title ('Original Signal Spectrum');

% Apply a 4kHz bandpass filter
cutoff_frequency = 4000;
audio_spectrum(f >= cutoff_frequency | f <= -cutoff_frequency) = 0;
filtered_signal = ifft(ifftshift(audio_spectrum));

% Plot the filtered signal in the frequency domain
audio_spectrum_filtered = fftshift (fft (filtered_signal));
f_filtered = fs_audio/2 * linspace(-1, 1, length (filtered_signal) ) ;
subplot (2,2,2);
plot (f_filtered, abs (audio_spectrum_filtered) / length (filtered_signal) );
title ('M(f) ');

% Calculate the time vector for the filtered signal
start_time = 0;
end_time = start_time + length (filtered_signal) / fs_audio;
time_filtered = linspace(start_time, end_time, length (filtered_signal))';
subplot (2,2,3);
plot (time_filtered, filtered_signal);
title('m(t)');

% Play the filtered signal
sound(abs(filtered_signal), fs_audio);
pause (2) ;

% Calculate modulation parameters
bandwidth_message = cutoff_frequency;
Fc = 100000;
amplitude_message = max(filtered_signal);
amplitude_carrier = 2 * amplitude_message;

% Resample the signal at 5 times the carrier frequency
filtered_signal_resampled = resample(filtered_signal, 5 * Fc, fs_audio) ;
fs_resampled = 5 * Fc;

% SSB-SC Modulation (Double Sideband Suppressed Carrier)
% Time vector for resampled signal
time_resampled = linspace(0, length (filtered_signal_resampled) / fs_resampled, length(filtered_signal_resampled) )'

% Generate DSB-SC modulated signal
carrier_wave = cos (2 * pi * Fc * time_resampled);
DSB_SC = filtered_signal_resampled .* carrier_wave;

% Generate SSB by filtering the lower sideband (LSB)


SSB_LSB = DSB_SC;
f_ssb = fs_resampled/2 * linspace(-1, 1, length(SSB_LSB));
audio_spectrum_ssb = fftshift(fft(SSB_LSB));
audio_spectrum_ssb (f_ssb >= Fc | f_ssb <= -Fc) = 0;
subplot (2,2,4);
plot (f_ssb, abs (audio_spectrum_ssb) / length(SSB_LSB));
title('LSB Spectrum') ;

% Coherent Detection of SSB-SC using ideal filter
demodulated_signal = SSB_LSB .* cos (2 * pi * Fc * time_resampled) ;
demodulated_fft = fftshift(fft(demodulated_signal));
demodulated_fft (f_ssb >= bandwidth_message | f_ssb <=- bandwidth_message) = 0;
demodulated_time = ifft(ifftshift(demodulated_fft));

% Plot the demodulated signal in the time domain
figure;
subplot (2,2,1);
plot (time_resampled, demodulated_time) ;
title ('Demodulated Signal using Ideal Filter waveform') ;

% Plot the demodulated signal in the frequency domain
signal_length_demodulated = length(demodulated_time);
demodulated_fft = fftshift(fft(demodulated_time));
f_demodulated = fs_resampled / 2 * linspace(-1, 1, signal_length_demodulated) ;
subplot (2,2,2);
plot (f_demodulated, abs (demodulated_fft) / signal_length_demodulated);
title ('Demodulated Signal using Ideal Filter Spectrum');

% Resample and play the demodulated signal
demodulated_resampled = resample (abs (demodulated_time), fs_resampled / 5, fs_resampled);
sound (abs (demodulated_resampled), fs_resampled / 5);
pause (2) ;

% Coherent Detection of SSB-SC using Butterworth filter
demodulated_butter = SSB_LSB .* cos (2 * pi * Fc * time_resampled) ;
[b_coeff, a_coeff] = butter(4, bandwidth_message * 2 / fs_resampled) ;
demodulated_butter = filtfilt (b_coeff, a_coeff, demodulated_butter);

% Plot the demodulated signal using Butterworth filter in the time domain
subplot (2,2,3);
plot (time_resampled, demodulated_butter);
title ('Demodulated Signal using Butterworth Filter waveform') ;

% Plot the demodulated signal using Butterworth filter in the frequency domain
signal_length_butter = length (demodulated_butter) ;
demodulated_fft_butter = fftshift(fft(demodulated_butter));
f_butter = fs_resampled / 2 * linspace(-1, 1, signal_length_butter) ;
subplot (2,2,4);
plot (f_butter, abs (demodulated_fft_butter) / signal_length_butter);
title ('Demodulated Signal using Butterworth Filter Spectrum') ;

% Resample and play the demodulated signal using Butterworth filter
demodulated_butter_resampled = resample (abs (demodulated_butter), fs_resampled / 5, fs_resampled) ;
sound (abs (demodulated_butter_resampled), fs_resampled / 5) ;
pause (2) ;

% Coherent detection with 0dB SNR
noisy_ssb_0db = awgn (SSB_LSB, 0) ;

demodulated_noisy_0db=noisy_ssb_0db .* cos (2 * pi * Fc * time_resampled) ;
demodulated_noisy_fft_0db = fftshift(fft(demodulated_noisy_0db));
demodulated_noisy_fft_0db(f_ssb >= bandwidth_message | f_ssb <=- bandwidth_message) = 0;
demodulated_noisy_0db =ifft(ifftshift(demodulated_noisy_fft_0db));

% Plot the demodulated signal with 0dB SNR in time domain
figure;
subplot (2,1,1);
plot(time_resampled, demodulated_noisy_0db);
title('0dB SNR Demodulated Signal waveform');

% Plot the demodulated signal with 0dB SNR in frequency domain
signal_length_noisy_0db = length(demodulated_noisy_0db);
demodulated_noisy_fft_0db = fftshift(fft(demodulated_noisy_0db));
f_noisy_0db = fs_resampled / 2 * linspace(-1, 1, signal_length_noisy_0db);
subplot (2,1,2);
plot(f_noisy_0db, abs(demodulated_noisy_fft_0db) / signal_length_noisy_0db);
title('0dB SNR Demodulated Signal spectrum');

% Resample and play the demodulated signal with 0dB SNR
demodulated_noisy_resampled_0db = resample(abs(demodulated_noisy_0db), fs_resampled / 5, fs_resampled) ;
sound(abs(demodulated_noisy_resampled_0db), fs_resampled / 5);
pause (2);

% Coherent detection with 10dB SNR
noisy_ssb_10db = awgn(SSB_LSB, 10);
demodulated_noisy_10db = noisy_ssb_10db .* cos (2 * pi * Fc * time_resampled);
demodulated_noisy_fft_10db = fftshift(fft(demodulated_noisy_10db));
demodulated_noisy_fft_10db(f_ssb >= bandwidth_message | f_ssb <=- bandwidth_message) = 0;

demodulated_noisy_10db = ifft(ifftshift(demodulated_noisy_fft_10db));

% Plot the demodulated signal with 10dB SNR in time domain
figure;
subplot (2,1,1);
plot(time_resampled, demodulated_noisy_10db);
title ('10dB SNR Demodulated Signal waveform');

% Plot the demodulated signal with 10dB SNR in frequency domain
signal_length_noisy_10db =length(demodulated_noisy_10db);
demodulated_noisy_fft_10db = fftshift(fft(demodulated_noisy_10db));
f_noisy_10db = fs_resampled / 2 * linspace(-1, 1, signal_length_noisy_10db);
subplot (2,1,2);
plot(f_noisy_10db, abs (demodulated_noisy_fft_10db) / signal_length_noisy_10db);
title('10dB SNR Demodulated Signal spectrum');

% Resample and play the demodulated signal with 10dB SNR
demodulated_noisy_resampled_10db = resample(abs(demodulated_noisy_10db), fs_resampled / 5, fs_resampled) ;
sound(abs (demodulated_noisy_resampled_10db), fs_resampled / 5);
pause (2);

% Coherent detection with 30dB SNR
noisy_ssb_30db = awgn(SSB_LSB, 30);
demodulated_noisy_30db = noisy_ssb_30db .* cos (2 * pi * Fc * time_resampled);
demodulated_noisy_fft_30db = fftshift(fft(demodulated_noisy_30db));
demodulated_noisy_fft_30db (f_ssb >= bandwidth_message | f_ssb <=- bandwidth_message) = 0;
demodulated_noisy_30db =ifft(ifftshift(demodulated_noisy_fft_30db));

% Plot the demodulated signal with 30dB SNR in time domain
figure;
subplot(2,1,1);
plot(time_resampled, demodulated_noisy_30db);
title('30dB SNR Demodulated Signal in waveform');

% Plot the demodulated signal with 30dB SNR in frequency domain
signal_length_noisy_30db=length(demodulated_noisy_30db);
demodulated_noisy_fft_30db=fftshift(fft(demodulated_noisy_30db));
f_noisy_30db = fs_resampled / 2* linspace(-1, 1, signal_length_noisy_30db);
subplot(2,1,2);
plot(f_noisy_30db, abs(demodulated_noisy_fft_30db) / signal_length_noisy_30db);
title('30dB SNR Demodulated Signal spectrum');

% Resample and play the demodulated signal with 30dB SNR
demodulated_noisy_resampled_30db = resample(abs(demodulated_noisy_30db), fs_resampled / 5, fs_resampled);
sound(abs(demodulated_noisy_resampled_30db), fs_resampled / 5);
pause (2) ;

% SSB-TC Generation
SSB_TC = amplitude_carrier .* carrier_wave + SSB_LSB;
signal_length_tc =length(SSB_TC);
audio_spectrum_tc=fftshift(fft(SSB_TC));
f_tc = fs_resampled / 2 * linspace(-1, 1, signal_length_tc);

% Plot SSB-TC in the frequency domain
figure;
subplot(2,1,1);
plot(f_tc, abs(audio_spectrum_tc) / signal_length_tc);
title('SSB-TC spectrum');

% Envelope Detection for SSB-TC
envelope_SSB_TC = abs(hilbert(SSB_TC));
subplot(2,1,2);
plot(time_resampled, SSB_TC);
hold on;
plot(time_resampled, -envelope_SSB_TC, 'g-', time_resampled, envelope_SSB_TC, '-g', 'LineWidth', 1.5);
title('SSB-TC waveform with Envelope Detector');
hold off;
ylim([-3 3]);
xlim([3 4]);

% Resample and play the envelope-detected signal
envelope_SSB_TC_resampled = resample(envelope_SSB_TC, fs_resampled / 5, fs_resampled);
sound (abs (envelope_SSB_TC_resampled), fs_resampled / 5);
pause (2);