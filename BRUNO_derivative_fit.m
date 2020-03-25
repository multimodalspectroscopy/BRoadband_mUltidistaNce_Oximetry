function least_square = BRUNO_derivative_fit(param,order,slope_diff, extinction,wave_start,wave_end,attenuation_slope_function,wavelengths,boundaries, distance,distance_max)

start_indx = find(wavelengths== wave_start);
end_indx = find(wavelengths == wave_end);
mua = param(1).*extinction(:,4) + log(10).* (param(2).*extinction(:,2) + param(3).*extinction(:,3)); %mu_a = WF*mua + ln(10)*e*c(HHb) + ln(10)*e*c(HbO2)
mus = param(4)*(wavelengths'*0.001).^(- param(5))';
if strcmp(boundaries,'far') == 1
    slope_model = feval(attenuation_slope_function, distance_max, distance, mua, mus);
else 
    slope_model = feval(attenuation_slope_function, mua, mus, distance);
end
slope_model_diff = diff(slope_model, order);        
difference = slope_model_diff(start_indx:end_indx) - slope_diff(start_indx:end_indx);
least_square = sum(difference.^2);
