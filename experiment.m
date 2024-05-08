%% define parameters
mus_low  = 0.01;
mus_high  = 20;
mua_low  = 0.001;
mua_high  = 2;
Ns = 10;
Na = 10;
% mus = logspace(log10(mus_low), log10(mus_high), Ns);
% mua = logspace(log10(mua_low), log10(mua_high), Na);
mus = linspace(mus_low, mus_high, Ns);
mua = linspace(mua_low, mua_high, Na);

%% run experiment
for mu_s = mus
    for mu_a = mua
        try
            test3d(mu_s, mu_a);
        end
        try
            test2d(mu_s, mu_a);
        end
        close all;
    end
end

%% read results
mu_expected2d = zeros(Ns, Na);
mu_reg2d= zeros(Ns, Na);
full_mu2d= zeros(Ns, Na);
single_exp_mu2d= zeros(Ns, Na);
full2d_mu2d= zeros(Ns, Na);

mu_expected3d = zeros(Ns, Na);
mu_reg3d= zeros(Ns, Na);
full_mu3d= zeros(Ns, Na);
single_exp_mu3d= zeros(Ns, Na);
alg_mu3d= zeros(Ns, Na);
full_point_mu3d= zeros(Ns, Na);
single_point_mu3d= zeros(Ns, Na);
sd0_point_mu3d= zeros(Ns, Na);

for i = 1:size(mus')
    for k = 1:size(mua')
        folder_3d = "3D";
        folder_2d = "2D";
        file_name = "mua" + string(mua(k)) + "mus" + string(mus(i)) + ".csv";
        full_file_name_2d = fullfile(folder_2d, "Results", file_name);
        full_file_name_3d = fullfile(folder_3d, "Results", file_name);
        data_array_2d = csvread(full_file_name_2d);
        data_array_3d = csvread(full_file_name_3d);
        
        
        mu_expected2d(i, k) = data_array_2d(1);
        mu_reg2d(i, k) = data_array_2d(2);
        full_mu2d(i, k) = data_array_2d(3);
        single_exp_mu2d(i, k) = data_array_2d(4);
        full2d_mu2d(i, k) = data_array_2d(5);
        
        mu_expected3d(i, k) = data_array_3d(1);
        mu_reg3d(i, k) = data_array_3d(2);
        full_mu3d(i, k) = data_array_3d(3);
        single_exp_mu3d(i, k) = data_array_3d(4);
        alg_mu3d(i, k) = data_array_3d(5);
        full_point_mu3d(i, k) = data_array_3d(6);
        single_point_mu3d(i, k) = data_array_3d(7);
        sd0_point_mu3d(i, k) = data_array_3d(8);


    end
end

%% analyse results
figure;
surf(mus, mua, mu_expected2d)
hold on
surf(mus, mua, mu_reg2d)
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("mu reg 2d")

figure;
surf(mus, mua, mu_expected2d)
hold on
surf(mus, mua, full_mu2d)
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("full mu 2d")


figure;
surf(mus, mua, mu_expected2d)
hold on
surf(mus, mua, single_exp_mu2d)
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("single exp mu 2d")


figure;
surf(mus, mua, mu_expected2d)
hold on
surf(mus, mua, full2d_mu2d)
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("full2d mu 2d")

figure;
surf(mus, mua, mu_expected3d)
hold on
surf(mus, mua, mu_reg3d)
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("mu reg 3d")

figure;
surf(mus, mua, mu_expected3d)
hold on
surf(mus, mua, full_mu3d)
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("full mu 3d")


figure;
surf(mus, mua, mu_expected3d)
hold on
surf(mus, mua, single_exp_mu3d)
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("single exp mu 3d")



figure;
surf(mus, mua, mu_expected3d)
hold on
surf(mus, mua, real(alg_mu3d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("alg mu 3d")


figure;
surf(mus, mua, mu_expected3d)
hold on
surf(mus, mua, real(full_point_mu3d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("full point mu 3d")


figure;
surf(mus, mua, mu_expected3d)
hold on
surf(mus, mua, real(single_point_mu3d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("single point mu 3d")

figure;
surf(mus, mua, mu_expected3d)
hold on
surf(mus, mua, real(sd0_point_mu3d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("sd0 point mu 3d")



%% diff results
figure;
surf(mus, mua, 20*log(mu_reg2d./mu_expected2d))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
title("mu reg 2d diff")
zlim([-6,6])
view(45,30)



figure;
surf(mus, mua, 20*log(full_mu2d./mu_expected2d))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("full mu 2d diff")
view(45,30)


figure;
surf(mus, mua, real(20*log(single_exp_mu2d./mu_expected2d)))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("single exp mu 2d diff")
view(45,30)

figure;
surf(mus, mua, 20*log(full2d_mu2d./mu_expected2d))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("full2d mu 2d diff")
view(45,30)

figure;
surf(mus, mua, 20*log(mu_reg3d./mu_expected3d))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("mu reg 3d diff")
view(45,30)

figure;
surf(mus, mua, 20*log(full_mu3d./mu_expected3d))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("full mu 3d diff")
view(45,30)

figure;
surf(mus, mua, real(20*log(single_exp_mu3d./mu_expected3d)))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("single exp mu 3d diff")
view(45,30)

figure;
surf(mus, mua, 20*log(real(alg_mu3d)./mu_expected3d))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("alg mu 3d diff")
view(45,30)

figure;
surf(mus, mua, real(20*log(real(full_point_mu3d)./mu_expected3d)))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("full point mu 3d diff")
view(45,30)

figure;
surf(mus, mua, real(20*log(real(single_point_mu3d)./mu_expected3d)))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("single point mu 3d diff")
view(45,30)

figure;
surf(mus, mua, real(20*log(real(sd0_point_mu3d)./mu_expected3d)))
hold on
surf(mus, mua, zeros(Ns,Na), 'FaceColor', 'r', 'FaceAlpha', 0.5); % Red plane with transparency
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
zlabel("dB")
zlim([-6,6])
title("sd0 point mu 3d diff")
view(45,30)


%% diff results
figure;
surf(mus, mua, 20*log(mu_reg2d./mu_expected2d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("mu reg 2d diff")

figure;
surf(mus, mua, 20*log(full_mu2d./mu_expected2d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("full mu 2d diff")


figure;
surf(mus, mua, real(20*log(single_exp_mu2d./mu_expected2d)))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("single exp mu 2d diff")


figure;
surf(mus, mua, 20*log(full2d_mu2d./mu_expected2d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("full2d mu 2d diff")

figure;
surf(mus, mua, 20*log(mu_reg3d./mu_expected3d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("mu reg 3d diff")

figure;
surf(mus, mua, 20*log(full_mu3d./mu_expected3d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("full mu 3d diff")

figure;
surf(mus, mua, real(20*log(single_exp_mu3d./mu_expected3d)))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("single exp mu 3d diff")

figure;
surf(mus, mua, 20*log(real(alg_mu3d)./mu_expected3d))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("alg mu 3d diff")

figure;
surf(mus, mua, real(20*log(real(full_point_mu3d)./mu_expected3d)))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("full point mu 3d diff")

figure;
surf(mus, mua, real(20*log(real(single_point_mu3d)./mu_expected3d)))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("single point mu 3d diff")

figure;
surf(mus, mua, real(20*log(real(sd0_point_mu3d)./mu_expected3d)))
xlabel("mus [1/mm]")
ylabel("mua [1/mm]")
title("sd0 point mu 3d diff")


