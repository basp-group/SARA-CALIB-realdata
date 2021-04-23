function [A, At, Gw, S] = op_p_nufft_time_calib(p, N, Nn, No, Ns)

% Create the nonuniform gridding matrix and fft operators to be used for
% parallel processing
%
% in:
% p{:}[2] - nonuniformly distributed frequency location points for each
%           cell member which will be treated in parallel
% N[2]    - size of the reconstruction image
% Nn[2]   - size of the kernels (number of neighbors considered on each direction)
% No[2]   - oversampled fft from which to recover the non uniform fft via
%           kernel convolution
% Ns[2]   - fft shift
%
% out:
% A[@]          - function handle for direct operator
% At[@]         - function handle for adjoint operator
% G{:}[:][:]    - convolution kernel matrix (small) associated with each
%               patch in the fourier plane
% W{:}          - mask of the values that contribute to the convolution
% Gw[:][:]      - global convolution kernel matrix
%originally written by Onose (2016), modified by Dabbech (2019)


if numel(N) ==1    
    [A, At, Gw, S] = op_nufft_time_calib(p, N, Nn, No, Ns);
else
    fprintf("\nThis is valid for 1D nufft used in calib!\n")
    
end

end