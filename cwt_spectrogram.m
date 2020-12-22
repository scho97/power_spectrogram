function [powermat,wtcoefs,freq,time] = cwt_spectrogram(time, X, Fs, wname, verbose)
%% Function 'cwt_spectrogram'
% This function computes and plots a continuous wavelet spectrogram (i.e. a
% magnitude spectrogram) of an input data. The current version only
% supports the single-trial CWT spectrogram; in other words, this code
% works only for an input data of a single vector array as of now.

% USAGE
% Full Input : cwt_spectrogram(time, X, Fs, wname, verbose)
% Example    : cwt_spectrogram([], data, 1024, 'Morse', 'expand')
%            : cwt_spectrogram(time, data, 1024)

% INPUT
%    Variable     Data Type                 Description
% 1. time         [1 x N vector]          : time array
%                                           1) []
%                                           2) vector array
% 2. X            [1 x N vector]          : input data
% 3. Fs           [integer Z > 0]         : sampling frequency rate
% 4. wname        [string]                : type of analysis wavelet
%                                         : 'Morse', 'amor' (default), 'bump'
% 5. verbose      [string]                : display CWT spectrogram
%                                           1) 'default' : the minimum and maximum scales will be
%                                                          determine automatically by the function 'cwt'
%                                           2) 'expand'  : uses 'cwtfreqbounds' to deterime the min and max 
%                                                          wavelet bandpass freqencies for given signal length, 
%                                                          sampling frequency, and wavelet

% OUTPUT
%    Variable     Data Type                                            Description
% 1. powermat     [length(freq) x length(time) matrix]               : power of 'wtcoefs'
% 2. wtcoefs      [length(freq) x length(time) complex matrix]       : continuous wavelet transform of the input signal for specified scales and wavelet
% 3. freq         [1 x N vector]                                     : frequency vector of CWT

% REFERENCE
% [1] MATLAB Documentation. https://www.mathworks.com/help/wavelet/ref/cwt.html
% [2] MATLAB Documentaiton. https://www.mathworks.com/help/wavelet/ref/cwtfreqbounds.html
% [3] Ahmet Taspinar. (2018). A guide ofr using the Wavelet Transform in Machine Learning. Retrieved from 
%     http://ataspinar.com/2018/12/21/a-guide-for-using-the-wavelet-transform-in-machine-learning/ on Dec 11, 2020.

% Note: As indicated by MATLAB,the plot uses a logarithmic frequency axis (in powers of 10) since frequencies in the CWT are logarithmic.

% Written by SungJun Cho, December 2020
% Last modified Dec 14, 2020
% Copyright (c) 2020 SungJun Cho
% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
% License: https://creativecommons.org/licenses/by-nc-sa/4.0/
%% Input Parameters
if nargin < 5
    verbose = 'default';
end
if nargin < 4
    wname = 'amor'; % Morlet Wavelet
end
%% Preprocess Data
% Create Time Array
if isempty(time)
    time = linspace(1, length(X)/Fs, length(X));
end
% Set Important Parameters
% dt = 1/Fs;
% Or alternatively: dt = time(2) - time(1);
%% Calculate CWT Spectrogram
switch verbose
    case 'default'
        [wtcoefs,freq] = cwt(X,wname,Fs);
    case 'expand'
        [minf,maxf] = cwtfreqbounds(numel(X),Fs);
        fb = cwtfilterbank('SignalLength', numel(X), 'Wavelet', wname, 'SamplingFrequency', Fs,'FrequencyLimits',[minf,maxf]);
        [wtcoefs,freq] = cwt(X,'FilterBank',fb);
        %figure(); freqz(fb);
end
powermat = abs(wtcoefs); % Or alternatively: abs(wtcoefs).^2;
%                                            wtcoefs.*conj(wtcoefs);
%% Plot CWT Spectrogram
% You can use the function 'surface' as suggested in MATLAB document, 
% but note that it is slower than 'imagesc'.
switch verbose
    case 'default'
        figure();
        imagesc(time,freq,imgaussfilt(powermat,1));
        shading flat; colormap('jet');
        ylim([1 max(freq)]);
        ftemp = flip(freq);
        ax = gca;
        ax.YTick = round(ftemp(1:10:end),3);
        set(ax,'yscale','log');
        axis xy;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title('Wavelet Spectrogram / Scalogram');
        c = colorbar;
        ylabel(c,'Coefficient Magnitude');
    case 'expand'
        figure();
        fBin = 15;
        freq_log = logspace(log10(minf),log10(maxf),fBin);
        %freq_log = flip(freq_log);
        imagesc(time,freq,imgaussfilt(powermat,1));
        shading flat; colormap('jet');
        ylim([1 max(freq_log)]);
        ax = gca;
        set(ax,'yscale','log');
        %ax.YTick = flip(freq_log);
        ax.YTick = freq_log;
        axis xy;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title('Wavelet Spectrogram / Scalogram');
        c = colorbar;
        ylabel(c,'Coefficient Magnitude');
end

% Display spectrogram properties
switch verbose
    case 'default'
        display_summary(wname,verbose,freq);
    case 'expand'
        display_summary(wname,verbose,freq,fb);
end
end

%% Accessory Functions
function display_summary(wname,verbose,freq,fb)
    if nargin < 4
        fb = [];
    end
    disp('*****************SUMMARY*****************');
    disp('CWT Spectrogram Parameters:');
    disp(['    Type of Wavelet: ' wname]);
    switch verbose
        case 'default'
            disp(['    Frequency Limits: ' num2str(min(freq)) 'Hz-' num2str(max(freq)) 'Hz']);
            disp(['    Number of Voices Per Octave: ' num2str(10)]);
        case 'expand'
            disp(['    Frequency Limits: ' num2str(fb.FrequencyLimits(1)) 'Hz-' num2str(fb.FrequencyLimits(2)) 'Hz']);
            disp(['    Number of Voices Per Octave: ' num2str(fb.VoicesPerOctave)]);
    end
    disp('*****************************************');
end