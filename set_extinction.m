function [Extinction_int,Waves_all,Wavelengths,DPF_wavelength_dependence] = set_extinction
% function to take extinction spectra and return them in micromoles per mm,
% also loads DPF 
% Current setting works for CYRIL, change files if required

DPF_wavelength_dependence = DPF_Lambda_Dependency_740to915;

Waves_all = csvread('PIXIS 512F wavelengths.txt'); %CYRIL pixel wavelengths
Wavelengths = [704:1:911]'; %Interpolation wavelengths
Water_ext = water_extinction_coeff;
Extinction = Wray_tissue_specific_extinction_coefficient_650to1042;
coeff = zeros([length(Wavelengths) 3]);
water_coeff = zeros([length(Wavelengths) 2]);
for i = 1:length(Wavelengths)
coeff(i,:) =  Extinction(find(Wavelengths(i) == Extinction(:,1)),1:3);
water_coeff(i,:) = Water_ext(find(Wavelengths(i) ==Water_ext(:,1)),1:2);
end
%water_coeff(:,2) = smooth(water_coeff(:,2),13,'sgolay',3);
water_abs = [water_coeff(:,1) water_coeff(:,2).*log(10)]; %absorption coefficient of water
ext_coeff = [coeff(:,1) coeff(:,3).*0.0000001 coeff(:,2).*0.0000001 water_abs(:,2)*0.1]; %the specific extinction coefficient of deoxy and oxy haemoglobin in microMoles per mm and the absorption coefficient of water per mm 
Extinction_int = zeros(length(Wavelengths),4);
Extinction_int(:,1) = ext_coeff(:,1);
for i = 2:4
   Extinction_int(:,i) = smooth(ext_coeff(:,i),13,'sgolay',3);
end