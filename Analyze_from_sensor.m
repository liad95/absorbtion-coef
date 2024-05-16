%% Analyzing straight from sensor
cutoff = 30;
cutoff_end = 90;
x_arr = -dx*(Nx-1)/2:dx:dx*(Nx-1)/2;
H_from_sensor = squeeze(sensor_data_rs(11,11,:));
x_arr = X(1,:,1);
lm = fitlm(x_arr(cutoff:cutoff_end), log(H_from_sensor(cutoff:cutoff_end)));
intercept = lm.Coefficients.Estimate(1);
mu_reg = -lm.Coefficients.Estimate(2);
simple_reg_energy = exp(-mu_reg*x_arr+ intercept);
fprintf("mu from simple linear regression %d\n", mu_reg)
figure;
plot(x_arr(cutoff:cutoff_end), H_from_sensor(cutoff:cutoff_end), "DisplayName","Real Data")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_energy(cutoff:cutoff_end), "DisplayName","Linear Regression")
title("Absorbtion for Different Regressions")
xlabel("Depth [mm]")
legend()
title("From Sensor")


figure;
plot(x_arr, H_from_sensor, "DisplayName","data")
hold on
plot(x_arr, log(H_from_sensor), "DisplayName","log")
hold on
line([x_arr(1)+(cutoff-1)*dx x_arr(1)+(cutoff-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
hold on
line([x_arr(1)+(cutoff_end-1)*dx x_arr(1)+(cutoff_end-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
xlabel("Depth [mm]")
legend()
title(["Absorbtion in middle ", name])
