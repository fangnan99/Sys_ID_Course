function periodogram = mypdg(F)

size_F = size(F);
length_F = size_F(1);
pdg_length = floor(length_F/2) +1;
answer = zeros(pdg_length, 1);
for ind = 1:1:pdg_length
    answer(ind) = abs(F(ind))*abs(F(ind)) / length_F;
end
periodogram = answer;
end