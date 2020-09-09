function F = mydft(f)
ii = sqrt(-1);
size_f = size(f);
length = size_f(1);
answer = zeros(length, 1);
for index = 1:1:length
    tmp = 0;
    for ind = 1:1:length
        tmp = tmp + f(ind) * exp(-ii*2*pi*(index-1)/length*(ind-1));
    end
    answer(index) = tmp;
end
F=answer;
end
