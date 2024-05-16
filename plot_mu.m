%% full sensor

sigmas = [0, 0.3, 0.5, 1];
mu_ideal = 2.184;
mu_expected = flip([1.884	2.037	2.133	2.184]);
mu_linear_valoMC = [2.156, 2.040, nan, nan];
mu_linear_kWave = [2.629, 2.813, 2.740, 4.921];
figure;
subplot(2, 1, 1);

plot(sigmas,mu_expected, 'mo--', DisplayName="mu expected value")
hold on;
plot(sigmas,mu_linear_valoMC, 'bo--', DisplayName="mu linear ValoMC")
hold on;
plot(sigmas,mu_linear_kWave, 'ko--', DisplayName="mu linear K-Wave")
hold on
yline(mu_ideal, 'r--', 'LineWidth', 2, DisplayName="mu ideal")
legend;
xlabel("sigma")
title("Mu vs Sigma - Full Sensor")



%% small sensor
sigmas = [0, 0.3, 0.5, 1];
mu_ideal = 2.184;
mu_expected = flip([1.884	2.037	2.133	2.184]);
mu_linear_valoMC = [2.147, 1.951, nan, nan];
mu_linear_kWave = [2.142, 2.284, 2.252, 3.726];
subplot(2, 1, 2);
plot(sigmas,mu_expected, 'mo--', DisplayName="mu expected value")
hold on;
plot(sigmas,mu_linear_valoMC, 'bo--', DisplayName="mu linear ValoMC")
hold on;
plot(sigmas,mu_linear_kWave, 'ko--', DisplayName="mu linear K-Wave")
hold on
yline(mu_ideal, 'r--', 'LineWidth', 2, DisplayName="mu ideal")
legend;
xlabel("sigma")
title("Mu vs Sigma - Small Sensor")

%% db
mu_linear_valoMC_db = 20*log(mu_linear_valoMC./mu_ideal);
mu_linear_kWave_db = 20*log(mu_linear_kWave./mu_ideal);
mu_expected_db = 20*log(mu_expected./mu_ideal);
figure;
plot(sigmas,mu_expected_db, 'mo--', DisplayName="mu expected value")
hold on;
plot(sigmas,mu_linear_valoMC_db, 'bo--', DisplayName="mu linear ValoMC")
hold on;
plot(sigmas,mu_linear_kWave_db, 'ko--', DisplayName="mu linear K-Wave")
hold on
legend;
xlabel("sigma")
ylabel("dB")
title("Mu vs Sigma - Small Sensor")


