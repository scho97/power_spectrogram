function [normSpec, muSpec, stdSpec] = normalize_spectrogram(type, method, Spec_f, Spec_t, Spec, Fs, verbose)
%% Function 'normalize_spectrogram'
% This function normalizes the input spectrogram at each frequency level.

% USAGE
% Full Input : [normSpec, muSpec, stdSpec] = normalize_spectrogram(type, method, Spec_f, Spec_t, Spec, verbose)
% Example    : [normSpec,~,~] = normalize_spectrogram('FT','zscore', Spec_f, Spec_t, Spec)
%            : [normSpec,~,~] = normalize_spectrogram('WT','minmax_scale',freq,time,powermat')

% INPUT
%    Variable       Data Type                                      Description
% 1. type           [string]                                     : type of spectrogram
%                                                                  1) FT: Fourier transfrom based
%                                                                  2) WT: Wavelet transfrom based
% 2. method         [string]                                     : time array
%                                                                  1) minmax_scale
%                                                                  2) minmax_mean
%                                                                  3) zscore
%                                                                  4) log (do not recommend for 'WT' as frequency is already in a logarithmic scale)
% 3. Spec_f         [1 x N vector]                               : input data
% 4. Spec_t         [1 x N vector]                               : sampling frequency rate
% 5. Spec           [length(Spec_t) x length(Spec_f) matrix]     : the lower threshold of frequency of interest
% 6. Fs             [integer Z > 0]                              : sampling frequency
% 6. verbose        [{0, 1, 2}]                                  : type of the plot
%                                                                  1) 0 - do not plot
%                                                                  2) 1 - plot without smoothing (not recommended for 'FT'-'log')
%                                                                  3) 2 - plot with smoothing

% OUTPUT
%    Variable       Data Type                                       Description
% 1. normSpec       [length(Spec_t) x length(Spec_f) matrix]      : the normalized spectrogram
% 2. muSpec         [1 x length(Spec_f) vector]                   : mean values of the spectrogram at each frequency
% 3. stdSpec        [1 x length(Spec_f) matrix]                   : standard deviation values of the spectrogram at each frequency


% Written by SungJun Cho, December 2020
% Last modified Dec 08, 2020
% Copyright (c) 2020 SungJun Cho
% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
% License: https://creativecommons.org/licenses/by-nc-sa/4.0/
%% Input Parameters
if nargin < 7
    verbose = 0;
end
%% Initialize Empty Arrays
muSpec = zeros(1,length(Spec_f));
stdSpec = zeros(1,length(Spec_f));
normSpec = zeros(length(Spec_t), length(Spec_f));
%% Mean and Standard Deviation of Spectrogram at Each Frequency
for i = 1:length(Spec_f)
    signal = Spec(:,i);
    muSpec(1,i) = mean(signal);
    stdSpec(1,i) = std(signal);
end
%% Normalization
switch method
    case 'minmax_scale'
        labelY = 'MinMax-Scale';
        for i = 1:length(Spec_f)
            signal = Spec(:,i);
            signal = (signal - min(signal)) ./ (max(signal)-min(signal));
            normSpec(:,i) = signal;
        end
    case 'minmax_mean'
        labelY = 'MinMax-Mean';
        for i = 1:length(Spec_f)
            signal = Spec(:,i);
            signal = (signal - muSpec(1,i)) ./ (max(signal)-min(signal));
            normSpec(:,i) = signal;
        end
    case 'zscore'
        labelY = 'Z-score';
        for i = 1:length(Spec_f)
            signal = Spec(:,i);
            signal = (signal - muSpec(1,i)) ./ stdSpec(1,i);
            normSpec(:,i) = signal;
        end
    case 'log'
        normSpec = pow2db(Spec);
        labelY = 'dB';
end
%% Plot the Normalized Spectrogram
switch type
    case 'FT'
        if verbose ~= 0
            figure();
            colormap('jet');
            if verbose == 1
                imagesc(Spec_t, Spec_f, normSpec');
            elseif verbose == 2
                imagesc(Spec_t, Spec_f, imgaussfilt(normSpec', 2));
            end
            axis xy;
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            c = colorbar;
            ylabel(c, ['Normalized Power (' labelY ')']);
            axis tight;
        end
    case 'WT'
        if verbose ~= 0
            figure();
            colormap('jet');
            if verbose == 1
                imagesc(Spec_t, Spec_f, normSpec');
            elseif verbose == 2
                imagesc(Spec_t, Spec_f, imgaussfilt(normSpec', 2));
            end
            ylim([1 max(Spec_f)+1]);
            shading flat;
            ax = gca;
            [minf,maxf] = cwtfreqbounds(numel(Spec_t),Fs);
            if round(maxf) == round(max(Spec_f))
                fBin = 15;
                freq_log = logspace(log10(minf),log10(maxf),fBin);
                ax.YTick = freq_log;
            else
                ftemp = flip(Spec_f);
                ax.YTick = round(ftemp(1:10:end),3);
            end
            set(ax,'yscale','log');
            axis xy;
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            c = colorbar;
            ylabel(c, ['Normalized Power (' labelY ')']);
        end
end
end