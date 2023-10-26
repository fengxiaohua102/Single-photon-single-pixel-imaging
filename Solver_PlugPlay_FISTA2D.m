function [x,convergence] = Solver_PlugPlay_FISTA2D(A,At,b,x0,opt)
% The objective function: F(x) = 1/2 ||y - hx||^2 + lambda |x|
% Input:
%   A: the forward operator function, input is a 3d cube, output a 2d
%   matrix
%   At: backward operator function, input is a 2D image, output a 3d matrix
%   b: degraded image (2D matrix)
%   x0: initialization (3D matrix)
%   opt.lambda: weight constant for the regularization term
% Output:
%   x: output image stack (3D)
%
% Author: Xiaohua Feng

fprintf(' - Running FISTA with fixed step size method\n');
lambda = opt.lambda;
maxiter = opt.maxiter;
vis = opt.vis;
xk = x0;
yk = xk;
tk = 1;
L0 = opt.step;
POSCOND= opt.POScond;
% k-th (k=0) function, gradient, hessian
objk  = func(x0,b,A,lambda);
fprintf('%6s %9s %9s\n','iter','f','sparsity');
fprintf('%6i %9.2e %9.2e\n',0,objk,nnz(xk)/numel(xk));
switch (opt.denoiser)
    case 'BM3D'
        denoise = @denoise_BM3D;
    case 'DnCNN'
        [ net, modelSigma ] = get_CNN_model(lambda);
        denoise = @(noise_x,sigma_n) denoise_DnCNN(net,noise_x,sigma_n);
    case 'BM4D'
        denoise = @denoise_BM4D;
    case 'Med2D'
        denoise = @denoise_Med2D;
    case 'VBM3D'
        denoise = @denoise_VBM3D;
    case 'ProxTV'
        denoise = @denoise_ProxTV;
    case 'ProxTV3D'
        denoise = @denoise_ProxTV3D;
    case 'ProxTV_Med'
        denoise = @denoise_ProxTVMed;
    case 'mixed'
        denoise = @denoise_mixed;
    case 'Prox_l1'
        denoise = @denoise_l1;
    case 'Prox_l1_GPU'
        denoise = @denoise_l1_GPU;
    otherwise
        denoise = @denoise_ProxTV;
end

for i = 1:maxiter
    i
    x_old = xk;
    y_old = yk;
    t_old = tk;
%     L0 = L0.*0.995;
    yg = y_old - 1/L0*(At(A(y_old)-b));

    zk = denoise(reshape(yg,size(x0)),lambda/L0); % denoise(yg,lambda/L0);
    fx = func(x_old,b,A,lambda);
    fz = func(zk,b,A,lambda);
    if(opt.monotone)
        if(fz<=fx)
            xk = zk;
        else
            xk = x_old;
        end
    end
    if(POSCOND)
        xk(xk<0)=0; % positiveness constraint
    end
    if( isfield(opt, 'ax') )
        set(opt.ax, 'CData', xk); drawnow;
    end
    tk = (1/20+sqrt(1/2+4*t_old*t_old))/2;
    yk = xk + (t_old-1)/tk*(xk-x_old)+(t_old)/tk*(zk-xk);
    if(POSCOND)
        yk(yk<0)=0; % positiveness constraint
    end
    if vis > 0
        fprintf('%6i %9.2e %9.2e\n',i,func(xk,b,A,lambda),nnz(xk)/numel(xk));
    end

    convergence(i) = fx;
    if(i>=3 && abs(convergence(end)-2*convergence(end-1)+convergence(end-2))/convergence(end-1)<=opt.tol )
        disp('Solution stagnates, exit ...');
        break;
    end
end
x = xk;
end

function norm_val = normNdMatrix(x,n)
    norm_val_temp = x.^n;
    norm_val = sum(norm_val_temp(:));
end

function Fx = func(xk,b,A,lambda)
    e = b - reshape(A(xk),size(b));
    % Fx = 0.5*normNdMatrix(e,2) + lambda*norm_tv3d(xk);
    Fx = 0.5*normNdMatrix(e,2);
end

function [img_estimated] = denoise_BM3D(noisy,th)
    norm_const = max(noisy(:));
    if(size(noisy,3)>1)
        noisy_3D = reshape(noisy,size(noisy,1),size(noisy,2),[]);
%         noisy_3D = prox_tv3d(noisy_3D,th);
        for K=1:size(noisy_3D,3)
            [~,x] = BM3D(1, noisy_3D(:,:,K)./norm_const, th*255, 'lc', 0); % NLmeansfilter
            x(isnan(x)) = 0;   noisy_3D(:,:,K) = x .* norm_const;
        end
        img_estimated = reshape(noisy_3D,size(noisy)) ;
    else
        [~,x] = BM3D(1, noisy./norm_const, th*255, 'lc', 0); % NLmeansfilter
        img_estimated = x*norm_const;
    end
end

function [img_estimated] = denoise_mixed(noisy,th)
    norm_const = max(noisy(:));
    if(size(noisy,3)>1)
        noisy_3D = reshape(noisy,size(noisy,1),size(noisy,2),[]);
%         noisy_3D = prox_tv3d(noisy_3D,th);
        for K=1:size(noisy_3D,3)
            [~,x] = BM3D(1, noisy_3D(:,:,K)./norm_const, th*255, 'lc', 0); % NLmeansfilter
            x(isnan(x)) = 0;   noisy_3D(:,:,K) = x .* norm_const;
            noisy_3D(:,:,K) = medfilt2(noisy_3D(:,:,K),[3,3]);
        end
        img_estimated = reshape(noisy_3D,size(noisy)) ;
    else
        [~,x] = BM3D(1, noisy./norm_const, th*255, 'lc', 0); % NLmeansfilter
        x(isnan(x)) = 0;        noisy = x.*norm_const;
%         img_estimated = prox_tv(noisy, th/1500);
        img_estimated = medfilt2(noisy,[3,3]);
    end

end


function x_filt = denoise_VBM3D(noisy,sigma_hat)
    norm_const = max(noisy(:));
    noisy_3D = reshape(noisy,size(noisy,1),size(noisy,2),[]);
    [~,x] = VBM3D_quite(noisy_3D./norm_const, sigma_hat*255); %VBM3D_quite
    x(isnan(x)) = 0.0;
    x_filt = x.*norm_const;
    x_filt = reshape(x_filt,size(noisy));
end

function x_filt = denoise_BM4D(noisy,sigma_hat)
    norm_const = max(noisy(:));
    noisy_3D = reshape(noisy,size(noisy,1),size(noisy,2),[]);
    [x,~] = bm4d(noisy_3D./norm_const, 'Gauss',sigma_hat*255,'lc',1,0); %VBM3D_quite
    x(isnan(x)) = 0.0;
    x_filt = x.*norm_const;
    x_filt = reshape(x_filt,size(noisy));
end


function img = denoise_Med2D(noisy,sigma_hat)
	img = medfilt2(noisy,[round(sigma_hat)+1,round(sigma_hat)+1]);%[5,5]
    img = imgaussfilt(img, 0.3);
end

function img = denoise_ProxTV(noisy,sigma_hat)
    if(size(noisy,3)>1)
        noisy_3D = reshape(noisy,size(noisy,1),size(noisy,2),[]);
%         noisy_3D = prox_tv3d(noisy_3D,sigma_hat);
        for K = 1:size(noisy_3D,3)
             noisy_3D(:,:,K) = prox_tv(noisy_3D(:,:,K),sigma_hat); %
        end
        img = reshape(noisy_3D,size(noisy));
    else
        img = prox_tv(noisy,sigma_hat);
    end
end
function img = denoise_ProxTV3D(noisy,sigma_hat)
    if(size(noisy,3)>1)
        noisy_3D = reshape(noisy,size(noisy,1),size(noisy,2),[]);
        noisy_3D = prox_tv3d(noisy_3D,sigma_hat);
        img = reshape(noisy_3D,size(noisy));
    end
end

function img = denoise_ProxTVMed(noisy,sigma_hat)
    if(size(noisy,3)>1)
        noisy_3D = reshape(noisy,size(noisy,1),size(noisy,2),[]);
        img = prox_tv3d(noisy_3D,sigma_hat);
        img = reshape(img,size(noisy));
    else
        img = prox_tv(noisy,sigma_hat);
    end
    img = medfiltn(img,3);
end

function img_estimated = denoise_l1(noisy,sigma_hat)
    img_estimated = prox_l1(noisy, sigma_hat);
end

function img_estimated = denoise_l1_GPU(noisy,sigma_hat)
    img_estimated = prox_l1_GPU(noisy, sigma_hat);
end

function img = medfiltn(noisy, nx)
    if(size(noisy,3)>1)
        noisy_3D = reshape(noisy,size(noisy,1),size(noisy,2),[]);
        for K = 1:size(noisy_3D,3)
%             noisy_3D(:,:,K) = imgaussfilt(noisy_3D(:,:,K),0.5); %
            noisy_3D(:,:,K) = medfilt2(noisy_3D(:,:,K),[nx,nx]); %
        end
        img = reshape(noisy_3D,size(noisy));
    else
        img = medfilt2(noisy,[nx,nx]);
    end
end

function img = denoise_DnCNN(net, noisy,sigma_hat)
    if(size(noisy,3)>1)
        noisy_3D = gpuArray(reshape(noisy,size(noisy,1),size(noisy,2),[]));
        img = gpuArray(zeros(size(noisy_3D)));
        for K = 1:size(noise_3D, 3)
            res = vl_simplenn(net, noisy_3D(:,:,K),[],[],'conserveMemory', true, 'mode', 'test');
            img(:,:,K) = noisy - res(end).x;
        end
        img = gather(reshape(img,size(noisy)));
    else
        res = vl_simplenn(net, gpuArray(noisy),[],[],'conserveMemory', true, 'mode', 'test');
        img = gather(gpuArray(noisy) - res(end).x ./50);
    end
end

function [ net, modelSigma ] = get_CNN_model( noiseSigma)
    folderModel = 'D:\MyAPPs\MatlabPackages\DnCNN\model';

    % load [specific] Gaussian denoising model
    modelSigma  = min(75,max(10,round(noiseSigma/5)*5)); %%% model noise level
    load(fullfile(folderModel,'specifics',['sigma=',num2str(modelSigma,'%02d'),'.mat']));
    %%%
    net = vl_simplenn_tidy(net);
    %%% move to gpu
    net = vl_simplenn_move(net, 'gpu') ;
end
