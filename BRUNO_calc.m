 function [StO2, coefficients, residual,residual_norm, sum_residual,score] = BRUNO_calc(slope, extinction, attenuation_slope_function, wavelengths, plotting, boundaries, distance, distance_max)
%% Run BRUNO and return oxygenation 
%
% INPUTS:
% slope: slope of attenuation against distance, W x 1 array, where W is number of wavelengths
% extinction: the extinction matrix, W x 4, first column wavelength, second
%       HHb, third HbO2, fourth the absorption coeff of water. If the output of
%       "set_extinction" is used, the extinction units are gonna be set
%       appropriately for this code, otherwise, check units
% attenuation_slope_function: function, the theoretical model for the attenuation
%       slope, generated using "generate_model"
% wavelengths: W x 1 array of wavelengths
% plotting: a double, either 1 or 0. Set to 1 if you want to plot the final
%       model agaisnt the data, set to 0 if not.
% boundaries: the boundary conditions for the minimisation procedure. 3 x
%       5 array, where the first row are the start values, second row are the lower
%       boundaries and third row are the upper boundaries. First column is water
%       fraction, second is HHb, third is HbO2, third and fourth are the
%       scattering coefficients from the exponential model; a and b resp., from mus =
%       a*lambda^(-b) where lambda is in micrometers.
% distance: a double, the light source-detector separation. If I use the
%       slope model where I only use one distance, this will be it, if I use the
%       minimal a maximal distance, set this as the minimal distance.
% distance_max: a double, only used if I use the slope model with two
%       different light source-detector separations. In that case, set this to
%       the maximal separation. Leave empty otherwise!
%
% OUTPUTS:
% StO2: the StO2 at the time point in %
% coefficients: the found model coefficients, first water fraction, then
%       HHb, HbO2, a and b 
% residual: 1 x (W - 1) array of squared residuals of the fit
% residual_norm: 1 x (W - 1) array of squared residuals of the normalised
%       fit
% sum_residual: 1 x 2 array, first number is the sum of residuals, second
%       is the sum of normalised residuals 
% score: a double, the score of the fit and the data
%


if nargin == 9
    bc = 'far';
else 
    bc = 'close';
    distance_max = 0;
end 

start= boundaries(1,:);  
LB = boundaries(2,:);  
UB =boundaries(3,:);  

slope_1stdiff = diff(smooth(slope,5)); %differentiate with step size 1nm

% Set wavelength range the fitting is performed at
wave_start = 710;
wave_end = 900;
options = optimset('Display','off','MaxIter',200000,'MaxFunEvals',200000,'TolX',1e-10,'TolFun',1e-10); %setting fitting options

% Performs the minimisation
Objective_function = @(param) BRUNO_derivative_fit(param,1,slope_1stdiff, extinction,wave_start,wave_end,attenuation_slope_function,wavelengths,bc,distance,distance_max);
[coefficients] = fminsearchbnd(Objective_function,start,LB,UB,options);

% Calculates mua and mus from the recovered coefficients for the plot
% and residual calculation
mua = coefficients(1).*extinction(:,4) + log(10).* (coefficients(2).*extinction(:,2) + coefficients(3).*extinction(:,3)); %mu_a = WF*mua + ln(10)*e*c(HHb) + ln(10)*e*c(HbO2)
mus = coefficients(4)*(wavelengths'*0.001).^(- coefficients(5))';
if strcmp(bc,'far') == 1 %evaluating my model
    model = feval(attenuation_slope_function, distance_max, distance, mua, mus);
else 
    model = feval(attenuation_slope_function, mua, mus, distance);
end
model_1stdiff = diff(model); 

% Model plot
if plotting == 1
    BRUNO_plot(model_1stdiff,slope_1stdiff,wavelengths,wave_start,wave_end) 
end

% Calculate StO2
StO2 = coefficients(3)/(coefficients(2)+coefficients(3)) * 100;

% Calculate residuals
residual(1,:) = (model_1stdiff- slope_1stdiff).^2;
residual_norm(1,:) = (model_1stdiff./(max(model_1stdiff))- slope_1stdiff./(max(model_1stdiff))).^2;
sum_residual(1) = sum(residual);
sum_residual(2) = sum(residual_norm);

% Calculate residuals only across HHb and water peak and use it to
% calculate the score
index_HHb = [find(wavelengths == 750):find(wavelengths == 770)];
index_water = [find(wavelengths == 825):find(wavelengths == 840)];
sum_hhb_residuals = sum(residual_norm(index_HHb));
sum_water_residuals = sum(residual_norm(index_water));
range =max(model_1stdiff(7:197)./(max(model_1stdiff(7:197)))) - min(model_1stdiff(7:197)./(max(model_1stdiff(7:197)))); %penalise the fit if it's too narrow 
score = sum_hhb_residuals.*sum_water_residuals*1./range;

