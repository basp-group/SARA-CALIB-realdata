function x = so_fft_adj_time_calib(X,N,No,scale)
x = prod(No).* (ifft(X)); 
x = x(1:N).*conj(scale);
end
