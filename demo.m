%% Load Data
% NOTE: Due to the current institutional restrictions, 'y_maze1_data000.mat' and 'y_maze1_data001.mat'
% files cannot be supplied to general public. Here, you should change the code and load your own data.
data1 = load('y_maze1_data000.mat');
data2 = load('y_maze1_data001.mat');
dataPFC1 = data1.EEG.data(1,:);
dataPFC2 = data2.EEG.data(1,:);
dataPFC = cat(2,dataPFC1,dataPFC2);
Fs = data1.EEG.srate; % sampling rate (Hz)
Nyq = Fs/2; % Nyquist frequency
%% Preprocess Data
% Apply Notch Filter: Eliminate 60 Hz Line Noise and its harmonics
filttype = 'iirnotch';

switch filttype
    case 'iirnotch'
        for f = 60:60:Nyq
            w = f/Nyq;
            Q = 30; % Q factor
            bw = w/Q;
            [b,a] = iirnotch(w,bw);
            dataPFC = filtfilt(b,a,dataPFC);
        end
    case 'butter_loword'
        for f = 60:60:Nyq
            [b,a] = butter(2,[f-1,f+1]/Nyq,'stop');
            dataPFC = filtfilt(b,a,dataPFC);
        end
    case 'butter_highord'
        for f = 60:60:Nyq
            [z,p,k] = butter(10,[f-1,f+1]/Nyq,'stop');
            [sos,g] = zp2sos(z,p,k);
            dataPFC = filtfilt(sos,g,dataPFC);
        end
end
% Detrend: Remove polynomial trend
dataPFC = detrend(dataPFC);
%% FFT-Hanning Spectrogram
lf = 0;
hf = 200;
window_size = 1024;
dt = 0.01;
[Spec_f, Spec_t, Spec] = hann_spectrogram([], dataPFC, Fs, lf, hf, window_size, dt);

method = 'log';
verbose = 2;
[normSpec, muSpec, stdSpec] = normalize_spectrogram('FT', method, Spec_f, Spec_t, Spec, Fs, verbose);
%% Multitaper Spectrogram
nw = 3;
ntp = 5;
M = 1024;
low_foi = 0;
high_foi = 200;
[Spec_t, Spec_f, Spec_taper_mean] = mtp_spectrogram([], dataPFC, Fs, nw, ntp, M, low_foi, high_foi);

method = 'minmax_scale';
verbose = 2;
[~, ~, ~] = normalize_spectrogram('FT', method, Spec_f, Spec_t, Spec_taper_mean, Fs, verbose);
%% Continuous Wavelet Transform (CWT) Spectrogram / Scalogram
wname = 'amor';
verbose = 'expand';
[powermat,wtcoefs,freq,time] = cwt_spectrogram([], dataPFC, Fs, wname, verbose);

method = 'minmax_scale';
verbose = 2;
Spec = powermat';
[~, ~, ~] = normalize_spectrogram('WT', method, freq, time, Spec, Fs, verbose);