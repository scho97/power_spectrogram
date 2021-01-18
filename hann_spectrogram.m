function [Spec_f, Spec_t, Spec] = hann_spectrogram(time, X, Fs, lf, hf, window_size, dt)
%% Function 'hann_spectrogram'
% This function computes and plots a FFT-Hanning spectrogram of an input
% data. The current version only supports the single-trial multitaper
% spectrogram (i.e. works only for an input data of a signle vector array
% as of now).

% USAGE
% Full Input : hann_spectrogram(time, X, Fs, lf, hf, window_size, dt)
% Example    : hann_spectrogram([], data, 1024, 0, 100, 1024, 0.1)

% INPUT
%    Variable       Data Type                    Description
% 1. time           [1 x N vector]             : time array
%                                                1) []
%                                                2) vector array
% 2. X              [1 x N vector]             : input data
% 3. Fs             [integer Z > 0]            : sampling frequency rate
% 4. lf             [integer Z >= 0 & =< Nyq]  : the lower threshold of frequency of interest
% 5. hf             [integer Z >= 0 & =< Nyq]  : the upper threshold of frequency of interest
% 6. window_size    [integer Z > 0]            : the length of Hanning window
%                                                default) 2^(nextpow2(Fs))
% 7. dt             [rational number Q > 0]    : time step by which the window slides (seconds)
%                                                default) 0.1

% OUTPUT
%    Variable         Data Type                                          Description
% 1. Spec_t           [1 x floor(diff(t_fft)/dt)+1 vector]               : vector of times for the spectrogram
% 2. Spec_f           [1 x length(lf:hf) vector]                         : vector of frequency bins for the spectrogram
% 3. Spec_taper_mean  [floor(diff(t_fft)/dt)+1 x length(lf:hf) matrix]   : matrix of calculated spectral power

% Written by SungJun Cho, December 2020
% Last modified Dec 07, 2020
% Copyright (c) 2020 SungJun Cho
% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
% License: https://creativecommons.org/licenses/by-nc-sa/4.0/
%% Input Parameters
if nargin < 7
    dt = 0.1;
end
if nargin < 6
    window_size = 2^(nextpow2(Fs));
end

% Create Time Array
if isempty(time)
    time = linspace(0, length(X)/Fs, length(X));
end

%% Calculate FFT-Hanning Spectrogram
t_fft = [time(1)+(((window_size*0.5)+1)/Fs) time(end)-(((window_size*0.5)+1)/Fs)];
t_vec = linspace(t_fft(1), t_fft(2), floor(diff(t_fft)/dt)+1);
num_time = length(t_vec);

% Getting sliding window indices
idx_collection = zeros(num_time,2);
for tIdx = 1:length(t_vec)
    idx_collection(tIdx,1) = find(time < t_vec(tIdx), 1, 'last') - window_size*0.5;
    idx_collection(tIdx,2) = find(time < t_vec(tIdx), 1, 'last') + window_size*0.5-1;
end

idx_diff = idx_collection(:,2) - idx_collection(:,1);
if idx_diff ~= (window_size-1)
    error('Error: Index computation incorrect!');
end

% Define frquency resolution
nfft = max(2^(nextpow2(window_size)), window_size);
df = Fs/nfft;
f = 0:df:Fs/2;

select_freq = find(f >= lf & f <= hf);
num_freq = length(select_freq);
Spec = zeros(num_time, num_freq);
Spec_f = f(select_freq);
Spec_t = time(floor(mean(idx_collection,2)));

hanning = hann(window_size)';

for time_idx = 1:length(t_vec)
    x_hann = hanning .* X(idx_collection(time_idx,1):idx_collection(time_idx,2));
    x_fft = fft(x_hann) / length(x_hann); % or alternatively: x_fft = fft(d,nfft) / Fs;
    Spec(time_idx,:) = abs(x_fft(select_freq));
    % or alternatively: abs_val = conj(x_fft).*x_fft; Spec(time_idx,:) = abs_val(select_freq,:);
end

%% Plot Multitaper Spectrogram
figure();
colormap('jet');
imagesc(Spec_t, Spec_f, imgaussfilt(Spec', 1));
%imagesc(Spec_t, Spec_f, imgaussfilt(pow2db(Spec'), 1));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
c = colorbar;
ylabel(c, 'Power ($\mu$V/Hz)', 'Interpreter', 'latex','FontSize',12);
axis tight;
title('Hanning Spectrogram')

% Display spectrogram properties
display_summary(df, dt, window_size, Fs);
end

%% Accessory Functions
function display_summary(df, dt, window_size, Fs)
    disp('*****************SUMMARY*****************');
    disp('FFT-Hanning Spectrogram Parameters:');
    disp(['    Frequency Bin: ' num2str(df) 'Hz']);
    disp(['    Window Step (Time Bin): ' num2str(dt) 's']);
    disp(['    Window Length: ' num2str(window_size/Fs) 's']);
    disp('*****************************************');
end
