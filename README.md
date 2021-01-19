[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC&20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Release](https://img.shields.io/github/release/scho97/power_spectrogram.svg)](https://github.com/scho97/power_spectrogram/releases/latest)

# Power Spectrogram
This repository features MATLAB scripts for computing and visualizing Fourier and wavelet transform-based power spectrograms. FT-based methods include `hann_spectrogram.m` and `mtp_spectrogram.m`, which produces a single taper power spectrogram using a Hanning window function and a multitaper power spectrogram that employs discrete prolate spheroidal sequences (DPSS), respectively. `cwt_spectrogram.m` is a WT-based method and performs continuous wavelet transform with a Morlet wavelet in default. Once the spectrogram is attained using one of these three methods, it can be normalized by frequency using `normalize_spectrogram.m`. 
Examples of signal preprocessing and power spectrogram generations are demonstrated in `demo.m`, which has to be supplied with the user's own data for the current version v1.1.

### Notes
The Python version is in progress and will be updated soon. If you are unfamiliar with signal processing and the codes appear unclear to you, please let me know. If there are any errors you noticed or any contents you would like me to add or remove, I would love to hear more about them. I can be reached by [email](mailto:scho20@uchciago.edu).

## Helpful Links
1. For the elimination of pink noise from the time-frequency data, refer to this [page](https://www.mathworks.com/matlabcentral/fileexchange/82470-pink-1-f-noise-elimination-for-time-frequency-data?s_tid=srchtitle).
2. For the tutorial on the multitaper spectrogram and another version of the code, visit [Prerau Lab](https://prerau.bwh.harvard.edu/multitaper/).

## License
This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/).  
Copyright (c) 2020 SungJun Cho [![CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
