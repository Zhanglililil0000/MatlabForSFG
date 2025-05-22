clear;
close all;

%% 导入实验数据
import_filepath = '\';
import_filename = 'WaterSSP.xlsx';
import_file = [import_filepath,import_filename];
RawData = readmatrix(import_file);

import_wavenumber = RawData(:,1);
import_intensity = RawData(:,2);


%% 基本实验参数
vis_angle_degree = 45; % 可见入射角
ir_angle_degree = 55; % 红外入射角
vis_wavelength = 532.1; % 可见波长

%% 中间参数
% 波长
import_wavelength = 1e7 ./ import_wavenumber;
sfg_wavelength = 1./(1./vis_wavelength + 1./import_wavelength);

% 入射反射角的弧度
vis_angle_rad = deg2rad(vis_angle_degree);
ir_angle_rad = deg2rad(ir_angle_degree);
sfg_angle_rad = asin(sfg_wavelength .* (sin(vis_angle_rad)./vis_wavelength +  sin(ir_angle_rad)./import_wavelength));

% 折射率
n_air = 1; % 空气折射率
n_quartz_vis = calculate_quartz_refractive_index(vis_wavelength); % 可见光折射率
n_quartz_ir = calculate_quartz_refractive_index(import_wavelength); % 红外折射率
n_quartz_sfg = calculate_quartz_refractive_index(sfg_wavelength); % SFG折射率


% 折射角
vis_ref_angle_rad = calculate_refraction_angle(vis_angle_rad, n_air, n_quartz_vis); % 可见光折射角
ir_ref_angle_rad = calculate_refraction_angle(ir_angle_rad, n_air, n_quartz_ir); % 红外折射角
sfg_ref_angle_rad = calculate_refraction_angle(sfg_angle_rad, n_air, n_quartz_sfg); % SFG折射角

% 相干长度(Unit: nm)
sfg_term = sqrt(n_quartz_sfg.^2 - sin(sfg_angle_rad).^2) ./ sfg_wavelength;
vis_term = sqrt(n_quartz_vis^2 - sin(vis_angle_rad)^2) / vis_wavelength;
ir_term = sqrt(n_quartz_ir.^2 - sin(ir_angle_rad).^2) ./ import_wavelength;
            
coherence_length = 1 ./(2 .* pi .* (sfg_term + vis_term + ir_term));  % 单位为nm

% 石英的菲涅耳因子
% SFG
lxx_sfg = quartz_fresnel(n_air, n_quartz_sfg, sfg_angle_rad, sfg_ref_angle_rad, 'xx');
lyy_sfg = quartz_fresnel(n_air, n_quartz_sfg, sfg_angle_rad, sfg_ref_angle_rad, 'yy');
          
% VIS
lxx_vis = quartz_fresnel(n_air, n_quartz_vis, vis_angle_rad, vis_ref_angle_rad, 'xx');
lyy_vis = quartz_fresnel(n_air, n_quartz_vis, vis_angle_rad, vis_ref_angle_rad, 'yy');

% IR
lxx_ir = quartz_fresnel(n_air, n_quartz_ir, ir_angle_rad, ir_ref_angle_rad, 'xx');
lyy_ir = quartz_fresnel(n_air, n_quartz_ir, ir_angle_rad, ir_ref_angle_rad, 'yy');


%% 计算石英的二阶极化率 单位(m^2/V)
% SSP
chi2_quartz_ssp = cos(ir_angle_rad) .* lyy_sfg .* lyy_vis .* lxx_ir .* coherence_length .* 1e-9 .* 1.6e-12;
chi2_quartz_ssp_Intensity = chi2_quartz_ssp.^2;

% PPP
chi2_quartz_ppp = cos(sfg_angle_rad) .* cos(vis_angle_rad) .* cos(ir_angle_rad) .* lxx_sfg .* lxx_vis .* lxx_ir .* coherence_length .* 1.6e-21;
chi2_quartz_ppp_Intensity = chi2_quartz_ppp.^2;

% SPS
chi2_quartz_sps = cos(vis_angle_rad) .* lyy_sfg .* lxx_vis .* lyy_ir .* coherence_length .* 1.6e-21;
chi2_quartz_sps_Intensity = chi2_quartz_sps.^2;

% PSS
chi2_quartz_pss = cos(ir_angle_rad) .* lxx_sfg .* lyy_vis .* lyy_ir .* coherence_length .* 1.6e-21;
chi2_quartz_pss_Intensity = chi2_quartz_pss.^2;


%% 得到chi2_eff_ssp
chi2_eff_ssp = import_intensity .* chi2_quartz_ssp_Intensity;

figure;
subplot(3,1,1);
plot(import_wavenumber,chi2_quartz_ssp_Intensity,"black");
ylabel('Quartz Chi_{eff} ^2');
xlabel('wavenumber');
xlim([1450 3900]);
subplot(3,1,2);
yyaxis left;
plot(import_wavenumber,import_intensity,"blue");
ylabel('Normalized SFG');
xlabel('wavenumber');
ylim([-1e-3 5e-3]);
xlim([1450 3900]);
hold on;
yyaxis right;
plot(import_wavenumber,chi2_eff_ssp,"red");
ylabel('Chi_{eff} ^2');
ylim([-1.2e-43 6e-43]);


%% 计算菲涅尔
SFGLyy = fresnel(n_air, sfg_wavelength, sfg_angle_rad, sfg_ref_angle_rad, 'yy');
VISLyy = fresnel(n_air, vis_wavelength, vis_angle_rad, vis_ref_angle_rad, 'yy');
IRLzz = fresnel(n_air, import_wavelength, ir_angle_rad, ir_ref_angle_rad, 'zz');

%% 得到Chi_yyz
YYZ_local_factor = (SFGLyy .*VISLyy .*IRLzz .*sin(ir_angle_rad)).^2;
chi_yyz = chi2_eff_ssp ./YYZ_local_factor;

subplot(3,1,3);
yyaxis left;
plot(import_wavenumber,import_intensity,"blue");
ylabel('Normalized SFG');
xlabel('wavenumber');
ylim([-1e-3 5e-3]);
xlim([1450 3900]);
hold on;
yyaxis right;
plot(import_wavenumber,chi_yyz,"red");
ylabel('Chi_{yyz} ^2');


figure;
subplot(2,1,1);
n_bulk_IR = bulk_refractive_index_interp(import_wavelength); % 水的折射率
plot(import_wavenumber,n_bulk_IR);
subplot(2,1,2);
plot(import_wavenumber,YYZ_local_factor);


% 结果导出
writematrix([import_wavenumber import_intensity chi2_eff_ssp chi_yyz],'watersspeffyyz.csv');
writematrix([import_wavenumber chi2_quartz_ssp_Intensity SFGLyy VISLyy YYZ_local_factor],'waterEffYYZfactors.csv')


%% 各项计算函数
% 计算折射率的函数
function n = calculate_quartz_refractive_index(wavelength)
    wavelength_um = wavelength ./ 1000;
    n_squared = 1.28604141 + 1.07044083 .* wavelength_um.^2 ./ (wavelength_um.^2 - 0.0100585997) + 1.10202242 .* wavelength_um.^2 ./ (wavelength_um.^2 - 100);
    n = sqrt(n_squared);
end


% 计算折射角的函数
function gamma = calculate_refraction_angle(incident_angle, n1, n2)
    gamma = asin(n1 .* sin(incident_angle) ./ n2);
end

% 计算石英菲涅耳因子
function quartzfresnel_result = quartz_fresnel(n1, n2, theta1, theta2, polarization)
 
    cos_theta1 = cos(theta1);
    cos_theta2 = cos(theta2);
        
    if polarization == 'xx'
        % Lxx = (2 * cosθ2) / (cosθ2 + n2 * cosθ1)
        numerator = 2 .* cos_theta2;
        denominator = cos_theta2 + n2 .* cos_theta1;
        quartzfresnel_result =  numerator ./ denominator;
    elseif polarization == 'yy'
        % Lyy = (2 * cosθ1) / (cosθ1 + n2 * cosθ2)
        numerator = 2 .* cos_theta1;
        denominator = cos_theta1 + n2 .* cos_theta2;
        quartzfresnel_result = numerator ./ denominator;
    end
end

% 差值得到体相的折射率
function n_interp = bulk_refractive_index_interp(wavelength)
    %WATER_REFRACTIVE_INDEX_INTERP 根据波数插值计算水的复折射率
    %   输入:
    %       wavenumber - 波数(cm-1),可以是标量或向量
    %   输出:
    %       n_interp - 折射率实部
    %       k_interp - 折射率虚部
    
    % 读取水的折射率数据
    data_path = 'RefractiveIndex\';
    ref_data = readmatrix([data_path, 'nH2O.csv']);
    
    % 原始数据列: 波长(μm), n
    lambda = ref_data(:,1);  % 波长(um)
    lambdanm = lambda .* 1000; % 波长(nm)
    n_data = ref_data(:,2);  % 折射率实部
    
    % 创建插值函数(使用样条插值保证平滑性)
    n_interp = interp1(lambdanm, n_data, wavelength, 'spline');
end


% 界面菲涅耳因子计算
function fresnel_result = fresnel(n1, wavelength, beta, gamma, polarization)
    
    n2 = bulk_refractive_index_interp(wavelength);
    % HXH Thesis
    % n_prime = sqrt((n1.^2 + n2.^2 + 4)./(2 .* (n1.^-2 + n2.^-2 + 1))); 
    % Consistency paper
    n_prime = sqrt((n2.^2 .*( n2.^2 + 5))./(4 .* n2.^2 + 2));

    if polarization == 'xx'
        % L_xx = ( 2 * n1 * cos(gamma)) / (n1 * cos(gamma) + n2 * cos(beta))
        numerator = 2 .* n1 .* cos(gamma);
        denominator = n1 .* cos(gamma) + n2 .* cos(beta);
        fresnel_result =  numerator ./ denominator;
    elseif polarization == 'yy'
        % L_yy = ( 2 * n1 * cos(beta)) / (n1 * cos(beta) + n2 * cos(gamma))
        numerator = 2 .* n1 .* cos(beta);
        denominator = n1 .* cos(beta) + n2 .* cos(gamma);
        fresnel_result =  numerator ./ denominator;
    elseif polarization == 'zz'
        % L_xx = ( 2 * n2 * cos(beta)) / (n1 * cos(gamma) + n2 *
        % cos(beta)) * (n1 / n') ^ 2
        numerator = 2 .* n2 .* cos(beta);
        denominator = n1 .* cos(gamma) + n2 .* cos(beta);
        fresnel_result =  numerator ./ denominator .* (n1 ./ n_prime) .^2;
    end

end


