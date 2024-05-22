radii = [10,20,30,40,49];
for radius=radii
    fprintf('run for %d', radius)
    KwaveSmallSensor_Analyze(radius);
    close all;
    clc;
end