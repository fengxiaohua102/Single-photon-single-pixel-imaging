function img_recon = SPC_T_ACC(meas_data, code, opt)
    % implement the adjoint operator of SPC
    % meas_data:[N_meas,1]
    % img_recon: [m,n]
    % code:[m,n,N_meas]

    [m,n,N_meas] = size(code);   % number of measurment
    meas_data = reshape(meas_data, [N_meas,1]);

    im_sum = zeros(m,n);
    for K_m = 1:N_meas
        % 1: apply SPC operator
        im_coded = meas_data(K_m) .* squeeze(code(:,:,K_m));
        im_sum = im_sum + im_coded;
    end
    % 2: apply transform basis operator
    img_recon = im_sum;
end
