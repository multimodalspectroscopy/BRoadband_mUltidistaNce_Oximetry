function BRUNO_plot(model_1stdiff,slope_1stdiff,wavelengths,wave_start,wave_end)
%% Function to plot the results of the fitting in one figure


start_indx = find(wavelengths == wave_start);
end_indx = find(wavelengths == wave_end-1);

figure
hold on
plot(wavelengths(start_indx:end_indx),model_1stdiff(start_indx:end_indx))
plot(wavelengths(start_indx:end_indx),slope_1stdiff(start_indx:end_indx))
legend ('Model','Data')
title('First derivative fit of attenuation slope')
xlabel('Wavelength [nm]')
ylabel('\partial^2 A/\partial\rho\partial\lambda  [mm^{-1}nm^{-1}]')
xlim([710 900])
box on