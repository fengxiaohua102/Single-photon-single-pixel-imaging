function F_res = conv_model_F(psf_in, x, opt)
% The convolutional forward model
% 1 psf_in: the normalized point spread function or simply the kernel
% 2 x: the input image/scene

Nx = opt.Nx; Ny = opt.Ny;
%x = reshape(x,[Nx,Ny]);
psf = zeros(size(x));
psf( (Nx-size(psf_in,1))/2+1:(Nx+size(psf_in,1))/2, (Ny-size(psf_in,2))/2+1:(Ny+size(psf_in,2))/2 ) = psf_in;
Hs = fft2( ifftshift(psf) );  % Compute 2D spectrum of the PSF

f_conv = real( ifft2( Hs.*fft2(x) ) );   % the forward convolution in FFT
% For 1D and sparse 2D detectors,simply extract the convolution results
switch(opt.mode)
    case '1D'
        F_res = f_conv(round(size(f_conv,1)/2),:);
    case 'Sparse_2D'
        N_stride = opt.stride;
        F_res = f_conv(1:N_stride:end,1:N_stride:end);
    otherwise
        F_res = f_conv;
end
end
