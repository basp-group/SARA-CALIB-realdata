function X = so_fft_time_calib(x,No,scale)
% oversampled FFT
X = fft(x.*scale,No);
end
