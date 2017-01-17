% Zeeman Experiment
% weighted fitting of a line to the energy vs magnetic field
% the raw data is the angle obtained from the video camera.
% You will supply an estimate of the angular uncertainty for each angle
% measurement. Some typical data has been supplied below in the variable
% x (array of magnetic field values, in Gauss)
% alpha1 (initial angle of an unshifted bright fringe)
% alpha (array of angles you measure as a function of the magnetic
% delta_alpha1 (uncertainly in the initial angle)
% delta_alpha (uncertainty in the shifted angles)
% JPS & MJM 2016
clear all; clf;
h=6.63e-34; c=3e8; lam=644e-9;
% enter the data.
x = [1287.18596104  1505.68749629  1658.80938095  1794.94411507 1948.12964917  2114.60223365  2233.34119528  2338.41449766]; %the indep.
%which we convert from Gauss to Tesla
x=x/10000;
alpha1=0.77; %angle of unshifted spectral line (degrees)
%convert to radians
alpha1=alpha1*pi/180;
%the next line is the angle of one of the shifted lines corresponding to
%the various magnetic fields
alpha = [0.818 0.828 0.839 0.845 0.85 0.861 0.871 0.877];
%the angle
%convert to radians
alpha=alpha*pi/180;
beta1=asin(sin(alpha1)/1.46); %this is the internal angle in the Fabry

beta=asin(sin(alpha)/1.46);
% now calculate the energies
del_lam_over_lam=(cos(beta)/cos(beta1))-1;
E = -(h*c/lam)*del_lam_over_lam;
%and convert to micro_eV
E=1e6*E/1.6e-19;
%Get the uncertainties in the energies. This follows from the
%uncertainties in the angles, alpha1 and alpha
delta_alpha1 = 0.005/2;
delta_alpha = 0.005/2.*ones([1,length(alpha)]); %angle uncertainty
%convert to radians;
delta_alpha1=delta_alpha1*pi/180; delta_alpha=delta_alpha*pi/180;
%and find the uncertainty in beta
delta_beta1=delta_alpha1*(cos(alpha1)/(1.46*cos(beta1)));
delta_beta=delta_alpha.*(cos(alpha)./(1.46*cos(beta)));

deltaE=-(h*c/lam)*sqrt(((sin(beta)/cos(beta1)).*delta_beta).^2 +...
 (cos(beta)*sin(beta1)./(cos(beta1)).^2.*delta_beta1).^2);
deltaE=1e6*deltaE/1.6e-19; %the uncertainty in the energy in micro_eV
%plot the data
errorbar(x,E,deltaE,'k.');
xlabel('Magnetic Field (T)');
ylabel('Energy Shift (\mu eV)');
hold on
%Now compute the coefficients for the weighted fit. See equations 8.37,38
%and 39 in Taylor
w=1./(deltaE.^2); %the weights
del=sum(w*sum(w.*(x.^2)))-(sum(w.*x))^2;
A=(sum(w.*(x.^2))*sum(w.*E)-sum(w.*x)*sum(w.*x.*E))/del;
B=(sum(w)*sum(w.*x.*E)-sum(w.*x)*sum(w.*E))/del;
Fit=A+B*x;
plot(x,Fit);
%What are the uncertainties in A and B?
%These are given by formulas quoted in problem 8.19
sig_A=sqrt(sum(w.*(x.^2))/del);
sig_B=sqrt(sum(w)/del);
disp(['the intercept is ' num2str(A) ' with uncertainty ' num2str(sig_A)])
disp(['the slope is ' num2str(B) ' with uncertainty ' num2str(sig_B)])