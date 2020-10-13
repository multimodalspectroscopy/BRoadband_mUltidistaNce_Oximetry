%% Analyse example data with BRUNO
% The example data is an attenuation spectrum collected in a homogeneous
% dynamic blood phantom during deoxygenation. 
% The data were collected at 30, 25, 20, 15 mm. 
% 
% The attenuation is scaled due to the way reference is collected.

% Requirements: MATLAB symbolic toolbox, MATLAB curve fitting toolbox and
% fminsearchbnd available at: https://uk.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon

%% Set constants 
% These are already provided: 
% - boundary conditions in the format [start upper lower] in the order water
%   fraction, HHb, HbO2 a and b 
%   a and b are from the scattering model mu_s' = a * wavelength/1000 ^ (-b)
% - wavelengths
% - extinction: wavelength, the extinction coefficient of HHb, HbO2 in microMoles per mm and water absorption
%   coefficient per mm

%% Calculate the slope of attenuation against distance, approximated with linear regression
matrix = repmat([ones(4,1), SD_separations'],1, size(atten_example',2)); warning('off','MATLAB:rankDeficientMatrix')
B = matrix\atten_example';
slope = B(2,:);

%% Load attenuation slope model, which will be later evaluated in BRUNO
model_slope = generate_model('ZBC', 0, 'attenuation_slope');

%% Run BRUNO 
% Returns StO2, the fit coefficients (WF, HHb, HbO2, a and b)
% See BRUNO_calc help for detailed description of inputs and outputs
[StO2, coefficients,~,~,~,~] = BRUNO_calc(slope, extinction, model_slope, wavelengths, 1, boundaries, mean([15 20 25 30]));

