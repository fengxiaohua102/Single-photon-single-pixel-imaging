%************
% This script performs the SPC image reconstruction into scattering medium
% using the dpeth from defocus idea
%************
clear all; clc; close all;
% pkg load image
m = 32; n = m; % spation resolution
mat_content = load([pwd '/ExpData/Letter5.mat']); %
meas_data = mat_content.meas_data;
meas_data = meas_data(1:2:end,:) - meas_data(2:2:end,:);
codes = hadamard(32*32);
codes = reshape(codes, m,n,[]);

figure; subplot(2,3,1); plot(meas_data(2,:)); title('Raw SPC data');

meas_data_ih = sum(meas_data,2);

im_ihada = hadamard(32*32)' * meas_data_ih;
im_ihada = reshape(im_ihada,m,n);

meas_data = sum(meas_data(:,100:220),2);
rot_angle = 90;

% preprocess the codes and measurement data
N_meas = size(codes,3);
codes_temp = reshape(codes, [m*n,N_meas]);
for K  = 1: size(codes_temp,1)
    codes_temp(K,:) = codes_temp(K,:) - mean(codes_temp(K,:) );
end
fac = 1.0/mean(vecnorm(codes_temp,2,2));
codes = reshape(codes_temp.*fac, [m,n,N_meas]);

% SPC solver setup
% 3: process the measurement
meas_data = meas_data - mean(meas_data(:));
meas_data = meas_data ./ max(meas_data(:));
meas_data_ih = meas_data_ih - mean(meas_data_ih(:));
meas_data_ih = meas_data_ih ./ max(meas_data_ih(:));

%% Reconstructin Step 1: SPC reconstruction
% 1: model setup
A = @(x) SPC_F_ACC(x, codes);
AT = @(x) SPC_T_ACC(x, codes);

% 2.reconstruction using the DAMP algorithm
max_egival = power_iter(A,AT,zeros(m,n))

% 3.FSITA recon
opt.tol = 1e-10;
opt.maxiter = 500;      % param for max iteration
opt.lambda = 1e-3;      % param for regularizing param
opt.vis = -1;
opt.denoiser = 'ProxTV'; % option of denoiser: BM3D, ProxTV,ProxTV_Med
opt.POScond = 1;         % positiveness contraint on the solution
global GLOBAL_useGPU;
GLOBAL_useGPU = 0;
opt.monotone = 1;
opt.step = 1.0*max_egival;  % step size

[recon_now,convergence] = Solver_PlugPlay_FISTA2D(A,AT,meas_data,zeros(m,n),opt);

% reconstruction using DC (time-integrated intensity) light 
[recon_DC,convergence] = Solver_PlugPlay_FISTA2D(A,AT,meas_data_ih,zeros(m,n),opt);

subplot(2,3,2); imagesc(fliplr(imrotate( recon_DC(2:end-2,2:end-2), rot_angle, 'crop' )) );
colormap(hot); title('DC recon.'); axis square off;
subplot(2,3,3); imagesc(fliplr(imrotate( recon_now(2:end-2,2:end-2), rot_angle , 'crop') ));
colormap(hot); title('SPC recon.'); axis square off;

%% Reconstructin Step 2: a: generate focal stack
blur_array = 2:0.2:5.0; 
Ndepth = length(blur_array);

spsf = zeros(m,n,Ndepth);

im_focal_stack = zeros(m,n,Ndepth);

% 1: Solver setup
opt_model.Nx = size(codes,1); 
opt_model.Ny = opt_model.Nx;
opt_model.mode = 'full2D';

opt.maxiter = 1000;         % param for max iteration
opt.lambda = 2e-3;          % param for regularizing param 5e-2 6e0*max_egival
opt.denoiser = 'BM3D';     % option of denoiser: BM3D, ProxTV,ProxTV_Med
opt.POScond = 1;            % positiveness contraint on the solution
opt.step = 1.0*max_egival;  % step size

for K_s = 1: Ndepth
    psf_deconv = fspecial("gaussian",m,blur_array(K_s));
    spsf(:,:,K_s) = norm1(psf_deconv);
%     psf_deconv = norm1(spsf(:,:,K_s));
    B = @(x) conv_model_F(psf_deconv, x, opt_model);
    BT = @(x) conv_model_T(psf_deconv, x, opt_model);
    disp('calculating step size...');
    max_egival = power_iter(B,BT,zeros(m,n))
    opt.step = 1.0*max_egival;  % step size
    [im_recon,convergence] = Solver_PlugPlay_FISTA2D(B, BT, recon_now, zeros(m,n), opt);

    im_focal_stack(:,:,K_s) = norm1(im_recon);
end

%% Reconstructin Step 2: b: focal stack filtering
im_focal_stack_flt = zeros(size(im_focal_stack)); BD_SIZE = 0;
im_focal_stack_flt(BD_SIZE+1:end-BD_SIZE, BD_SIZE+1:end-BD_SIZE,: ) = im_focal_stack(BD_SIZE+1:end-BD_SIZE, BD_SIZE+1:end-BD_SIZE,: );
[~,im_focal_stack_flt] = VBM3D(im_focal_stack_flt,15);
im_focal_stack_flt(isnan(im_focal_stack_flt)) = 0;
Image_FM = zeros(1,length(blur_array));

% measureing the image quality using contrast evaluation
for K = 1:size(im_focal_stack_flt,3)
    Image_FM(K) = sum(sum(FocusMeasure(imgaussfilt(im_focal_stack_flt(:,:,K),2.0), 'SML')));

end
Image_FM = imgaussfilt(medfilt1(Image_FM,3),3);

[~,im_idx] = max(Image_FM);
fprintf('\nIdentified index for refocusing is %d \n',im_idx);

subplot(2,3,5); plot(Image_FM,'LineWidth',2);  axis square; title('FocalMeasure')
subplot(2,3,6); imagesc(fliplr((imrotate( im_focal_stack_flt(:,:,im_idx), rot_angle, 'bicubic', 'crop' )))); title('Deconv. recon.');
axis square off;
