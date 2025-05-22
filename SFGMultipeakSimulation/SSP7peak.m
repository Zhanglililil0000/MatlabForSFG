close all;
clear;

%% for Multipeak SFG Spectra Simulation

%% Peak Parameter
NR_SSP_real = -0.0651*1e-20;
NR_SSP_imag = 0;

Omegas = [1665.2 3209.5 3298.3 3445 3558.0 3665.0 3700.50];
gammas = [45.4 51.1 45 188 98 50.6 14.03];
SSPAmps = [1.49 -1.33 -0.54 -31.9 10.5 2.8 2.42];
SSPAmps = SSPAmps.*1e-20;

%% Spectra Calculation
IR_range = linspace(1500,3900,10000);
ZeroBaseLine = zeros(10000,1);
PeakNum = length(SSPAmps);
SSPAmp = NR_SSP_real + NR_SSP_imag .* 1i;


for q = linspace(1,PeakNum,PeakNum)
    SSPAmp = SSPAmp + LorAmp(SSPAmps(q), Omegas(q), gammas(q), IR_range);
end

SSP_intensity = (abs(SSPAmp)).^2;
SSP_real = real(SSPAmp);
SSP_imag = imag(SSPAmp);


%% Result Plot
figure;
set(gcf, 'Position', [100, 100, 600, 800]);

% Intensity
subplot(2,1,1);
plot(IR_range,SSP_intensity);
title("SSP Intensity");
xlabel("wavenumber");


% Real and Imaginary Part
subplot(2,1,2);
hold on;
plot(IR_range,SSP_real,"red");
plot(IR_range,SSP_imag,"blue");
plot(IR_range,ZeroBaseLine,"black--");
title("SSP Real and Imaginary Part");
xlabel("wavenumber");
legend("Real","Imaginary",'location','northwest');


% %% plot different component peak
% figure;
% hold on;
% set(gcf, 'Position', [100, 200, 900, 600]);
% for q = linspace(1,PeakNum,PeakNum)
%     SSPPeaks(:,q) = LorAmp(SSPAmps(q), Omegas(q), gammas(q), IR_range);
%     plot(IR_range,imag(SSPPeaks(:,q)),"--");
% end
% plot(IR_range,SSP_imag,"black");
% 
% figure;
% hold on;
% set(gcf, 'Position', [100, 200, 900, 600]);
% for q = linspace(1,PeakNum,PeakNum)
%     SSPPeaks(:,q) = LorAmp(SSPAmps(q), Omegas(q), gammas(q), IR_range);
%     plot(IR_range,abs(SSPPeaks(:,q)).^2,"--");
% end
% plot(IR_range,SSP_intensity,"black");


% %% Export Simulation Result
% writematrix([IR_range',SSP_real',SSP_imag'],"SSPRealImag.csv");
% writematrix([IR_range',PPP_real',PPP_imag'],"PPPRealImag.csv");
% writematrix([IR_range',SSP_intensity'],"SSPIntensity.csv");
% writematrix([IR_range',PPP_intensity'],"PPPIntensity.csv");
% writematrix([IR_range',imag(SSPPeaks)],"SSPImcomponent.csv");
% writematrix([IR_range',imag(PPPPeaks)],"PPPImcomponent.csv");
% writematrix([IR_range',real(SSPPeaks)],"SSPRecomponent.csv");
% writematrix([IR_range',real(PPPPeaks)],"PPPRecomponent.csv");
% writematrix([IR_range',abs(SSPPeaks).^2],"SSPIntcomponent.csv");
% writematrix([IR_range',abs(PPPPeaks).^2],"PPPIntcomponent.csv");

