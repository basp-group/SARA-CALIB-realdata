function X = so_fft2_shift(x, No, scale)
% Computes and oversampled and scaled FFT2 with the scale factors
% precomputed from nufft_init
%
% in:
% x[:][:]  - input image
% No[2]    - overscale size
% scale[:] - scale parameters precomputed by nufft_init
%
% out:
% X[:]     - FFT2 coefficients


% apply scaling factors
x = x .* real(scale); 

% oversampled FFT 
X = fftshift(fft2(x, No(1), No(2)));
X = X(:);

end