% By Andrew John Buggee

clear variables
%% ---- TWO STREAM MODEL ---

% This code assumes the following:
%   (1) a semi-infinite cloud, with a hard boundary at the top
%   and extending to infinity in the other direction. 
%   (2) We assume photons can only travel up or down. Therefore the two end
%   states are (a) the photon is absorbed in the cloud or (b) the photon
%   scatters back out the cloud top

% Bohren and Clothiaux list a solution to this infinite two-stream problem
% in their book, starting on page 264. The reflection out the cloud top,
% the irradiance that we measure, is defined on page 265, and called R.


% Compute the reflectance out of the cloud top across many wavelengths. We
% need the asymmetry parameter and the absorption coefficient

% compute the single scattering albedo and asymmetry parameter using mie
% theory for some wavelength w and some radius r. Also, specify whether to
% use a monodispersed distribution or a gamma droplet distribution

droplet_distribution = 'gamma';
effective_radius = 5:5:10;                      % microns - effective droplet radii
%wavelength = 300:20:2000;                   % nm - wavelength
wavelength = 350:10:2300;                   % nm - wavelength range of the HySICS spectrometer

mie_properties = zeros(length(wavelength),8,length(effective_radius));
legend_str = cell(1,length(effective_radius));


for rr = 1:length(effective_radius)
    mie_properties(:,:,rr) = interp_mie_computed_tables([wavelength',repmat(effective_radius(rr),length(wavelength),1)],droplet_distribution,false);
    
    legend_str{rr} = ['$r_{e}$ = ',num2str(effective_radius(rr)),' $\mu m$'];
end

single_scattering_albedo = reshape(mie_properties(:,6,:),length(wavelength),[]);
g = reshape(mie_properties(:,7,:), length(wavelength),[]);

%complex_index_refraction = mie_properties(:,5,:);

%absoprtion_coefficient = 4*pi*complex_index_refraction./(wavelength'./1e3);         % microns^(-1)

% ----- using the Towmey and Bohren approximation for absorption ---
% R = (sqrt(1 - g.*(1 - 0.85*effective_radius*absoprtion_coefficient)) - sqrt(0.85*effective_radius*absoprtion_coefficient))./...
%     (sqrt(1 - g.*(1 - 0.85*effective_radius*absoprtion_coefficient)) + sqrt(0.85*effective_radius*absoprtion_coefficient));


% ----- Using the mie estimated values for single scattering ----
R = (sqrt(1 - g.*single_scattering_albedo) - sqrt(1 - single_scattering_albedo))./...
    (sqrt(1 - g.*single_scattering_albedo) + sqrt(1 - single_scattering_albedo));


% plot the results of reflectivity versus wavelength for several
% different droplet sizes

figure; plot(wavelength, R)
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Reflectivity','Interpreter','latex')
title(['Two-Stream Reflectivity: semi-infinite ', droplet_distribution,'-dispersed cloud'],'Interpreter','latex')
grid on; grid minor
legend(legend_str, 'Interpreter','latex','Location','best','FontSize',18)
set(gcf,'Position',[0 0 1200, 800])


%% Least Squares Solution - Can I retrieve r?

% no, this didn't work. Why not?
% Something about the assumptions of the two equations I'm equalizing must
% be wrong. Check the assumptions of the Two-stream reflectivity I'm using,
% and compare them with those of the Twomey Bohren approximation for
% absorption.

index_wavelength = [2, 50, 80];         % wavelengths to use

absorption_coefficient = 4*pi*reshape(mie_properties(:,5,:), length(wavelength),[])./(wavelength'./1e3);  % microns^(-1)


G = [ones(length(index_wavelength),1), -0.85* absorption_coefficient(index_wavelength)];         % linear model

d = (G' * G)^(-1) *G' * R(index_wavelength);

%% How does the two stream multiple scattering solutions compare to Twomey and Bohren?

% We solved the two-stream multiple scattering reflectivity of a
% semi-infinite cloud. Twomey and Bohren do the same using H-functions to
% solve for the absorption due to normal incident flux. These should be similar.

radius2plot = 10;
index2plot = effective_radius==radius2plot;

solar_zenithAngle = 0;                  % deg

abs_TB = sqrt((1 - single_scattering_albedo(:,index2plot))./(1 - g(:,index2plot).*single_scattering_albedo(:,index2plot))) .* H(cosd(solar_zenithAngle), single_scattering_albedo(:,index2plot),0)';

% In a semi-infinte layer, only absorption and reflection out the cloud top
% exist as the two photon end states.

R_TB = 1 - abs_TB;

% plot the results of two stream reflectivity and the reflectivity
% estiamted using the absorption approximation by Twomey and Bohren

figure; plot(wavelength, R(:,index2plot))
hold on;
plot(wavelength, R_TB)
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Reflectivity','Interpreter','latex')
title('Reflectivity of a semi-infinite cloud','Interpreter','latex')
grid on; grid minor
legend('Two-Stream','Twomey-Bohren', 'Interpreter','latex','Location','best','FontSize',18)
set(gcf,'Position',[0 0 1200, 800])