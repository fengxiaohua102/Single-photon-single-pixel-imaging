function T_res = conv_model_T(psf_in, x, opt)
% The convolutional adjoint model
% 1 psf: the normalized point spread function or simply the kernel
% 2 x: the camera measurement

% For 1D and sparse 2D detectors,simply extract the convolution results
Nx = opt.Nx; Ny = opt.Ny; %N_stride = opt.stride;
psf = zeros(size(x));
psf( (Nx-size(psf_in,1))/2+1:(Nx+size(psf_in,1))/2, (Ny-size(psf_in,2))/2+1:(Ny+size(psf_in,2))/2 ) = psf_in;

x_img = zeros(Nx,Ny);

switch(opt.mode)
    case '1D'        
        x_img(round(Nx/2),:) = x;
    case 'Sparse_2D'
        x_img(1:N_stride:end,1:N_stride:end) = x;
    otherwise 
        x_img = x;   
end

Hs_conj = conj( fft2( ifftshift(psf) ));       % Compute 2D spectrum
T_res = real( ifft2( Hs_conj.*fft2(x_img) ) ); % The adjoint operator

end