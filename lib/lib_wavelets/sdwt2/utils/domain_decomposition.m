function rg = domain_decomposition(nchunks, N)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% q : index of the facet
% nchunks : number of factes
% N : size of the domain

% Adapted from my Julia codes
splits = round(linspace(0, N, nchunks + 1));
rg = zeros(nchunks, 2);

for q = 1:nchunks
   rg(q, :) = [(splits(q) + 1), splits(q + 1)]; 
end

end