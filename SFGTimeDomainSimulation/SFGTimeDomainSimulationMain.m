close all;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        electric field definition              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the time domine, unit: ps
tn = 2^15; 
t1 = linspace(1, tn, tn) / 1000;                     % 卷积之前的时间坐标
t2 = linspace(1, 2 * tn - 1, 2 * tn - 1) / 1000;     % response function和红外电场卷积之后的时间坐标

% define IR light field
A_IR = 1;                                            % amplitude of IR
Time_IR = 0.40;                                       % pulse with of IR
wn_IR = 2900;                                        % center wavenumber of IR, Unit cm-1

[EIR1,IRprofie1] = IRlight(A_IR, Time_IR, wn_IR, t1);
[EIR2,IRprofie2] = IRlight(A_IR, Time_IR, wn_IR, t2);

% define VIS light field
A_vis = 1;                                           % amplitude of visible
dt_vis = 0;                                          % center time position of visible
Time_vis = 90;                                       % pulse with of IR
wl_vis = 532;                                        % unit: nm

[Evis,visprofile] = VISlight(A_vis, dt_vis, Time_vis, wl_vis, t2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         vibration modes definition            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_NR = -0.062;                                           % Non-resonance term
Psi_NR = 0;
Aq = [1.93, 1.09];                               % Amplitude of mode q
Omegaq = [2916.88, 2919.66];                    % position of mode q
Gammaq = [3.61, 2.93];                            % Lorentzian width (half)
dOmegaq = [0, 0];                                % Gaussian-broadened factor

% response function
R_chi2 = ResponseFunction(A_NR, Psi_NR, Aq, Omegaq, Gammaq, dOmegaq,t1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  FID process                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P1, Re1, Im1, Isfg1] = SFGFID(R_chi2, EIR1, Evis);

% FID reocess
figure;
set(gcf, 'unit', 'centimeters', 'position', [10 18 20 8]);
plot(t2,Re1);
% xlim([0,100]);
hold on;
plot(t2,visprofile);
title("FID process");
legend("FID","Visible pulse profile");
xlabel("time/ps");

% FID reocess, intensity
figure;
set(gcf, 'unit', 'centimeters', 'position', [10 3 40 12]);
subplot(1,2,1);
plot(t2,Isfg1);
% xlim([0 100]);
ylim([0 1]);
hold on;
plot(t2,visprofile);
title("FID process");
legend("FID Intensity","Visible pulse profile");
xlabel("time/ps");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Frequency domain                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ISFG,IRfreqprofile, wn, f] = time2freq(P1, t2, tn, EIR2, wl_vis);

subplot(1,2,2);
plot(wn,ISFG ./ 1.5E4);
xlim([2750,3000]);
hold on;
plot(f,IRfreqprofile ./ 1E2);
xlim([2750,3000]);
title("SFG spectra");
legend("SFG Intensity","IR profile");
xlabel("wavenumber/cm-1");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Output                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Functions                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 定义红外电场
function [EIR,IRprofie] = IRlight(A_IR, Time_IR, wn_IR, t)
    c = 0.03;
    c_IR = Time_IR / (2 * (log(2))^0.5); % IR pulse width parameter
    OmegaIR = 2 .*pi .* c .* wn_IR; % convert to frequency, Unit 1/ps

    EIR = A_IR .* exp(-1.386 .* (t / c_IR).^2 ) .* exp(1i .* OmegaIR .* (t));
    IRprofie = A_IR .* exp(-1.386 .* (t / c_IR).^2 );
end


%% 定义可见电场
function [Evis,visprofile] = VISlight(A_vis, dt_vis, Time_vis, wl_vis, t2)
    c = 0.03;
    c_vis = Time_vis ./ (2 * (log(2)).^0.5);
    Omegavis = 2 .*pi ./ wl_vis .* c .* 1e7; % convert to frequency, Unit 1/ps

    visshape = exp(-1.386 .* ((t2 - dt_vis) / c_vis).^2 ); % visible时间上的形状 
    % visshape = exp(-1.386 .* ((t2 - dt_vis) / c_vis).^2 ) .* (1+0.2.*cos(0.8.*t2))./1.2; % visible时间上调制后的形状

    Evis = A_vis .* visshape .* exp(1i .* Omegavis .* ((t2 - dt_vis)));
    visprofile = A_vis .* visshape;
end


%% Response Function SFG响应函数
function R_chi2 = ResponseFunction(A_NR, Psi_NR, Aq, Omegaq, Gammaq, dOmegaq,t)
    
    ModeNum = length(Aq);                           % count the number of modes
    c = 0.03;                                       % light speed, cm/ps
    Freqq = 2 .*pi .* c .* Omegaq;                  % 将波数换算成频率
    T2q = 1 ./ (2 .* pi .* c .* Gammaq);            % using Gammaq to calculate T2q
    
    % Response Function R2
    R_chi2 = abs(A_NR) .* exp(1i .* Psi_NR) .* dirac(t); % dirac is delta function in Matlab
    for n = 1:ModeNum
        R_chi2 = R_chi2 - 1i .* heaviside(t) .* Aq(n) .* exp(-1i .* Freqq(n) .*t) .* exp(-t ./ T2q(n)) .* exp(-dOmegaq(n).^2 .* t .^ 2 ./ 2);
    end
    R_chi2 = R_chi2 ./ sum(Aq);
end


%% FID过程
function [P, ReP, ImP, Isfg] = SFGFID(R_chi2, EIR, Evis) 

    % interaction between molecule and IR, the electronic polarization
    P_IR = conv(R_chi2, EIR);
    P = P_IR .* Evis;

    % output
    ReP = real(P);
    ImP = imag(P);
    Isfg = abs(P) .^ 2;

end


%% 时域到频域的变换
function [ISFG, IRfreqprofile, wn, f] = time2freq(P, t, tn, EIR, wlvis)
    
    f = t * (1000 / (2 * tn * 0.03)) * 1000;
    wn = - f + 1 / wlvis * 1e7;

    SFG=abs(fft(P)) .^ 2;
    IRfreqprofile = abs(fft(EIR)) .^ 2;
    ISFG = SFG ./ IRfreqprofile; % normalization

end
