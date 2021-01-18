function [Spec_t, Spec_f, Spec_taper_mean] = mtp_spectrogram(time, X, Fs, nw, ntp, M, low_foi, high_foi, plot_psd, check_opt, plot_tp, segment_yn, t_start, t_end)
%% Function 'mtp_spectrogram'
% This function computes and plots a multitaper spectrogram of an input
% data. The current version only supports the single-trial multitaper
% spectrogram (i.e. works only for an input data of a single vector array 
% as of now).

% USAGE
% Full Input   : mtp_spectrogram(time, X, Fs, nw, ntp, M, low_foi, high_foi, plot_psd, check_opt, plot_tp, segment_yn, t_start, t_end)
% Example      : mtp_spectrogram([], data, 1024, 2.5, 4, 1024, 0, 100)

% INPUT
%    Variable     Data Type                   Description
% 1. time         [1 x N vector]            : time array
%                                             1) []
%                                             2) vector array
% 2. X            [1 x N vector]            : input data
% 3. Fs           [integer Z > 0]           : sampling frequency rate
% 4. nw           [integer Z > 0]           : the time-half bandwidth for the tapers
%                                           : 2.5, 3, 3.5, 4 recommended
% 5. ntp          [integer Z > 0]           : the number of tapers
% 6. M            [integer Z > 0]           : the window length of DPSS tapers
%                                           : Fs, Fs/2, Fs*2 recommended
% 7. low_foi      [integer Z > 0 & =< Nyq]  : the lower threshold of frequency of interest
% 8. high_foi     [integer Z > 0 & =< Nyq]  : the upper threshold of frequency of interest
% 9. plot_psd     [boolean / logical]       : display multitaper power spectral density estimates
%                                            1) true or 1
%                                            2) false or 0
% 10. check_opt    [string]                  : check optimal number of tapers based on user's selected time-half bandwidth
%                                            1) []
%                                            2) 'method1'
%                                            3) 'method2'
% 11. plot_tp      [boolean / logical]       : display DPSS tapers
%                                            1) true or 1
%                                            2) false or 0
% 12. segment_yn  [boolean / logical]       : segment input data
%                                            1) true or 1
%                                            2) false or 0
% 13. t_start     [integer Z > 0]           : time (s) in which the segment starts
% 14. t_end       [integer Z > 0]           : time (s) in which the segment ends

% OUTPUT
%    Variable         Data Type                                                      Description
% 1. Spec_t           [1 x floor(diff(t_mt)/dt)+1 vector]                          : vector of times for the spectrogram
% 2. Spec_f           [1 x length(low_foi:high_foi) vector]                        : vector of frequency bins for the spectrogram
% 3. Spec_taper_mean  [floor(diff(t_mt)/dt)+1 x length(low_foi:high_foi) matrix]   : matrix of calculated spectral power

% REFERENCE
% [1] Prerau M.J., Brown R.E., Bianchi M.T., Ellenbogen J.M., Purdon P.L. (2016). Sleep Neurophysiological Dynamics
%     Through the Lens of Multitaper Spectral Analysis. Physiology, 32(1): 60-92. DOI: 10.1152/physiol.00062.2015.
%     Copyright 2019 Michael J. Prerau, Ph.D., Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (adapted the idea of parallel computing)
% [2] Percival D.B., Walden A.T. (1993). Spectral Analysis for Physical Applications. UK: Cambridge University Press.
% [3] Chronux Analysis Software. http://chronux.org/chronuxFiles/Documentation/chronux/spectral_analysis/helper/dpsschk.html
% [4] Chronux Analysis Software. http://chronux.org/chronuxFiles/Documentation/chronux/spectral_analysis/continuous/mtspectrumc.html
% [5] MATLAB Documentation. https://www.mathworks.com/help/signal/ref/dpss.html

% Written by SungJun Cho, October 2020
% Last modified Dec 07, 2020
% Copyright (c) 2020 SungJun Cho
% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
% License: https://creativecommons.org/licenses/by-nc-sa/4.0/
%% Input Parameters
if nargin < 12
    segment_yn = false;
end
if nargin < 9
    plot_psd = true;
    check_opt = 'method1';
    plot_tp = true;
    segment_yn = false;
end

% Initialize Parallel Computation
if isempty(gcp) % Get current parallel pool
    parpool; % Create parallel pool on cluster
end
%% Preprocess Data
% Create Time Array
if isempty(time)
    time = linspace(0, length(X)/Fs, length(X));
end

% Segment the Time Series
if segment_yn
    X = X(t_start*Fs:t_end*Fs);
    time = linspace(t_start, t_end, (t_end-t_start)*Fs+1);
    
    figure();
    plot(time, X);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Origianl Signal');
end
%% Construct DPSS Tapers
dpss_params = {nw, ntp};

% Plot Multitaper Power Spectral Density Estimate
if plot_psd
    w = linspace(0,1,Fs);
    figure();
    pmtm(X, dpss_params, w*(Fs/2), Fs, 'DropLastTaper', false);
end

% Create DPSS Tapers
[dpss_win, eigvals] = dpss(M, dpss_params{1}, dpss_params{2});
[~, num_win] = size(dpss_win);

% Check for the optimal number of tapers
if ~isempty(check_opt)
    switch check_opt
        case 'method1'
            if floor(2*nw-1) < ntp
                opt = floor(2*nw-1);
                warning(['The optimal number of tapers for the chosen time-half bandwidth is: ', num2str(opt)]);
            end
        case 'method2'
            figure();
            stem(1:length(eigvals), eigvals, 'filled');
            hold on;
            plot(1:length(eigvals), 0.99*ones(length(eigvals),1));
            ylim([0 1.2]);
            opt = find(eigvals>0.99,1,'last');
            title('Proportion of Energy in [-w, w] of k-th Slepian Sequence');
            if ntp ~= opt
                warning(['The optimal number of tapers for the chosen time-half bandwidth is: ', num2str(opt)]);
            end
    end
end

dpss_win = dpss_win * sqrt(Fs); % for the normalization of the tapers (cf. [3])

% Plot DPSS tapers
if plot_tp
    figure();
    hold on;
    for i = 1:num_win
        txt = ['Win #', num2str(i)];
        plot(dpss_win(:,i), 'DisplayName', txt);
    end
    xlabel('Samples');
    xlim([0 M]);
    title(['DPSS, M=', num2str(M), ', NW=', num2str(nw)]);
    legend show;
end
%% Calculate The Multitaper Spectrogram
window_size = M; 

%%%%% Time Resolution (Bin) %%%%%
dt = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_mt = [time(1)+(((window_size*0.5)+1)/Fs) time(end)-(((window_size*0.5)+1)/Fs)];
t_vec = linspace(t_mt(1), t_mt(2), floor(diff(t_mt)/dt)+1);

% Getting sliding window indices
idx_collection = zeros(length(t_vec),2);
for tIdx = 1:length(t_vec)
    idx_collection(tIdx,1) = find(time < t_vec(tIdx), 1, 'last') - window_size*0.5;
    idx_collection(tIdx,2) = find(time < t_vec(tIdx), 1, 'last') + window_size*0.5-1;
end

idx_diff = idx_collection(:,2) - idx_collection(:,1);
if idx_diff ~= (window_size-1)
    error('Error: Index computation incorrect!');
end

%%%%%%%%%%%%%%%% Frquency Resolution (Bin) %%%%%%%%%%%%%%%%
nfft = max(2^(nextpow2(window_size)), window_size);
df = Fs/nfft;
f = 0:df:Fs/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

select_freq = find(f >= low_foi & f <= high_foi);
num_freq = length(select_freq);
Spec = zeros(length(t_vec),num_freq,ntp);
Spec_f = f(select_freq);
Spec_t = time(floor(mean(idx_collection,2)));

data_segments = zeros(length(t_vec), window_size);
for time_idx = 1:length(t_vec)
    data_segments(time_idx,:) = X(idx_collection(time_idx,1):idx_collection(time_idx,2));
end

parfor time_idx = 1:length(t_vec)
    X_seg = data_segments(time_idx,:);
    X_mt = repmat(X_seg', 1, ntp) .* dpss_win;
    x_mt = fft(X_mt) / length(X_mt); % or alternatively: x_mt = fft(X_mt,nfft) / Fs;
    Spec(time_idx,:,:) = abs(x_mt(select_freq,:));
    % or alternatively: abs_val = conj(x_mt).*x_mt; Spec(time_idx,:,:) = abs_val(select_freq,:);
end

Spec_taper_mean = mean(Spec, 3);

%% Plot Multitaper Spectrogram
figure();
colormap('jet');
imagesc(Spec_t, Spec_f, Spec_taper_mean');
%imagesc(Spec_t, Spec_f, pow2db(Spec_taper_mean'));
axis xy
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%yticks(low_foi:10:high_foi);
c = colorbar;
ylabel(c, 'Power ($\mu$V/Hz)', 'Interpreter', 'latex','FontSize',12);
axis tight;
title('Multitaper Spectrogram');

% Display spectrogram properties
display_summary(df, dt, M, Fs, nw, ntp);
end

%% Accessory Functions
function display_summary(df, dt, M, Fs, nw, ntp)
    disp('*****************SUMMARY*****************');
    disp('Multitaper Spectrogram Parameters:');
    disp(['    Spectral Resolution: ' num2str(2*nw/(M/Fs)) 'Hz']);
    disp(['    Frequency Bin: ' num2str(df) 'Hz']);
    disp(['    Window Step (Time Bin): ' num2str(dt) 's']);
    disp(['    Window Length: ' num2str(M/Fs) 's']);
    disp(['    Time Half-Bandwidth Product: ' num2str(nw)]);
    disp(['    Number of Tapers: ' num2str(ntp)]);
    disp('*****************************************');
    % The spectral resolution corresponds to the bandwidth of the main lobe in
    % the spectral estimate, which specifies the minimum distance between peaks
    % that can be resolved [1].
end
