mag_grid = -5:0.1:5;
otput_total = zeros(length(mag_grid),1);
for i = 1:length(mag_grid)
    u = mag_grid(i) * ones(4000,1);
    y = HS2019_SysID_midterm_p2_system_sim(19951482, u);
    otput_total(i) = mean(y(1001:4000));
end
plot(mag_grid, otput_total);

