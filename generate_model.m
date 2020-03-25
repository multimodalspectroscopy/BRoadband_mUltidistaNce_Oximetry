function Model = generate_model(BC, Sep_assum, Quantity)
% Function to generate a model to use in BF or BRUNO 
%
% INPUTS:
% BC -  boundary conditions: 'ZBC'(zero boundary) or 'EBC'(extrapolated boundary)
%       ZBC reflectance taken from Lindkvist(2013), Spectroscopy Letters 46:4
%       EBC reflectance taken from Kienle(1997),j. Opt. Soc. Am. A 14:1
%   
%       NOTE: both models are applied with the assumption that rho^2 >> z_0^2 and that mua << mus
%
% Sep_assum - do we assume that the separation between the detectors is
%       negligible compared to the separation between the detectors and the light
%       source?
%       If it is neglibigle: Sep_assum = 0
%       If it is not negligible: Sep_assum = 1
%
%       In the case of a larger separation, the slope of attenuation is
%       evaluated as the attenuation at the two distant light source
%       separations divided by the difference in the separations.
%       For more info, see Eq. 19 in Scholkmann(2014) Physiol. Meas. 35
%
% Quantity - Specify which model you would like:
%       'reflectance' OR 'attenuation' OR 'attenuation_slope'
% 
% OUTPUT:
% Model - the desired model saved in a matlab function

% EXAMPLE:
% If I want to get an EBC model of the attenuation slope at long
% separations:
% Model = generate_model('EBC', 0, 'attenuation_slope')
%

%% CHECK INPUT
if (strcmp(BC, 'EBC') || strcmp(BC, 'ZBC'))  ~= 1
    BC = input('Please specify boundary conditions. Type "EBC" for extrapolated boundary conditions or "ZBC" for zero boundary conditions. Submit with pressing Return.');
end

if ((Sep_assum == 1) || (Sep_assum == 0)) ~= 1
    Sep_assum = input('Please specify your assumptions on the detector separation. If the separation between detectors is negligible compared to the separation between the detector and the light source, type 0. Else, type 1. Submit with pressing Return.');
end

if (strcmp(Quantity, 'reflectance') || strcmp(Quantity, 'attenuation') || strcmp(Quantity, 'attenuation_slope')) ~= 1
    Quantity = input('Please specify the model you require. Type: "reflectance" OR "attenuation" OR "attenuation_slope". Submit with pressing Return.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(BC, 'ZBC') == 1
%% ZERO BOUNDARY CONDITIONS

    syms mus mua rho
    z0 = 1/mus;
    mueff = sqrt(3*mua*mus);
    Reflectance = z0 * mueff * exp(-mueff*rho)/(2*pi*rho^2);
    Reflectance_function = matlabFunction(Reflectance);
    Attenuation = 1/log(10)*(mueff*rho + 2*log(rho) - log((z0*mueff)/(2*pi)));
    Attenuation_function = matlabFunction(Attenuation);
    
    if Sep_assum == 0
        Attenuation_slope = 1/(log(10))*(sqrt(3*mus*mua) + 2/rho);
        Attenuation_slope_function = matlabFunction(Attenuation_slope);
        
    elseif Sep_assum == 1
        syms dl ds
        Attenuation_slope = 1/(log(10))*(sqrt(3*mus*mua) + 2*((log(dl/ds))/(dl-ds)));
        Attenuation_slope_function = matlabFunction(Attenuation_slope);
        
    end
    
else
%% EXTRAPOLATED BOUNDARY CONDITIONS

    syms x0(rho,mua,mus) x1(rho,mua,mus) 
    z0 = 1./mus;
    D = 1./(3*(mua+mus));
    zb = (1+0.493)./(1-0.493)*2.*D; % from Kienle(1997),j. Opt. Soc. Am. A 14:1. Valid for biological tissue.
    r1 = rho^2; % r1^2 = z0^2 = rho^2 but z0^2 is negligible compared to rho^2
    r2 = (z0 + 2*zb)^2 + rho^2;
    mueff = sqrt(3*mua*mus);
    Reflectance = 1./(4*pi)*(z0*(mueff + 1./sqrt(r1)) * (exp(-mueff.*sqrt(r1))./r1) + (z0 + 2*zb).*(mueff + 1/sqrt(r2)) .* (exp(-mueff.*sqrt(r2))./r2));
    Reflectance_function = matlabFunction(Reflectance);
    Attenuation = - log10(Reflectance); 
    Attenuation_function = matlabFunction(Attenuation);
        
    if Sep_assum == 0
        Attenuation_slope =  diff(Attenuation_function,rho); 
        Attenuation_slope_function = matlabFunction(Attenuation_slope);
        
    elseif Sep_assum == 1
        syms dl ds
        Attenuation_slope = (feval(Attenuation_function,mua,mus,dl) - feval(Attenuation_function,mua,mus,ds))/(dl-ds);
        Attenuation_slope_function = matlabFunction(Attenuation_slope);
        
    end
end

%% PASS RESULT

if strcmp(Quantity, 'reflectance')
    Model = Reflectance_function;
elseif strcmp(Quantity, 'attenuation')
    Model = Attenuation_function;
else 
    Model = Attenuation_slope_function;
end

end
