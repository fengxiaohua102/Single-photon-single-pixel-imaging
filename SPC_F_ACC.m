function meas_data = SPC_F_ACC(image_object, code)
    % This is the accelerated version
    % implement the forward operator for CLIP-SPC
    % image_object: [m,n]
    % code:[m,n,N_meas]
    % meas_data: [N_meas,1]

    [m,n,N_meas] = size(code);
    image_object = reshape(image_object,[m,n]);
    meas_data = zeros(N_meas,1);
    for K_m = 1:N_meas
        % 2: apply SPC operator
        im_coded = image_object .* squeeze(code(:,:,K_m));
        meas_data(K_m) = sum(im_coded(:));
    end
end
