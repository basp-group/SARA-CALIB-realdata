function rg = domain_decomposition_overlap(nchunks, N, d)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% q : index of the facet
% nchunks : number of factes
% N : size of the domain
% d : number of overlapping pixels

% Adapted from my Julia codes
splits = round(linspace(0, N, nchunks + 1));

r = rem(d, 2);
rg = zeros(nchunks, 2);
for q = 1:nchunks
    if q > 1
        id_start = (splits(q) + 1) - floor(d/2) - r;
    else
        id_start = (splits(q) + 1);
    end

    if q < nchunks
        id_end = splits(q + 1) + floor(d/2);
    else
        id_end = splits(q + 1);
    end

    rg(q, :) = [id_start, id_end];
end

end
