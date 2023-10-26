function [depth_map, Img_ALLFOCUS] = depth_from_focus(Img_Stack, Method)
% ******
% perform depth from focus 
% Img_Stack is the refocused image at different depths [Nx,Ny,Ndepth]
% thresh: the thresholding that mask out low intensity region
% NOTE: the Img_Stack shoud be aligned!
% ******

[Nx,Ny,N_STEP] = size(Img_Stack);
ImgFM_Stack = zeros(Nx,Ny,N_STEP);
% 1: compute the foucs measure at corresponding depths [Nx,Ny, Ndepth]
for K = 1:N_STEP
    img_flt =  imgaussfilt(Img_Stack(:,:,K),3.0); %
    im_fm = FocusMeasure(img_flt,Method);
%     im_fm = medfilt2(im_fm, [11,11]);
    im_fm = prox_tv(im_fm, 1e-3);
    ImgFM_Stack(:,:,K) = im_fm;
end

% 2: Select the best depth/slop based on the best focus-measure
% first form the all-in-focus image by selecting the best-in-focus stack index
[~, idx] = max(ImgFM_Stack,[],3);

idx = prox_tv(double(idx)./N_STEP, 1e-0);
idx = round(idx.*N_STEP);
% idx = round( medfilt2(idx,[9,9]));
idx(idx==0)=1;

Img_Org_Reshape = reshape(Img_Stack, Nx*Ny,N_STEP);
Img_ALLFOCUS = zeros(Nx*Ny,1);

for K = 1: size(Img_ALLFOCUS,1)
    Img_ALLFOCUS(K) = Img_Org_Reshape(K,idx(K));
end
Img_ALLFOCUS = reshape(Img_ALLFOCUS,Nx,Ny);
% mask = (Img_ALLFOCUS<thresh*max(Img_ALLFOCUS(:)));
depth_map = idx;
% depth_map(mask) = NaN;
end