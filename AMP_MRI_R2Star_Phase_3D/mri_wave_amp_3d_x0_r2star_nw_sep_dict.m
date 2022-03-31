    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 1: Setting up parameters %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ESPIRiT parameters
    ncalib = 24;        % size of central k-space, use 24 calibration lines to estimate senstivity map, make sure to do full sampling in the AC region

    % dataset parameters
    sx = 320;       % size along x-direction
    sy = 320;       % size along y-direction
    sz = 208;       % size along z-direction
    Nc = 32;        % the number of channels (coils)
    Ne = 4;         % the number of echoes
    echo_time = [7.32, 16, 24.7, 33.39].'; % echo time

    LN_num = 10;                        % computing thread
    LN = maxNumCompThreads( LN_num );   % set the largest number of computing threads

    % generate Poisson sampling pattern accroding to poisson disk
    sampling_rate = 0.1;        % the sampling rate
    sampling_rate_str = '20';
    sampling_pd_interval = 2;   % the minimum distance between any two locations, needs to be tuned for best performance
    num_samples = round(sy*sz*sampling_rate);   % the number of samples
    
    % for this dataset, sx direction is fully sampled, undersampling takes place in sy-sz plane
    % full_sampling_loc should be a mask of size sy by sz, the mask contains 0-1 values 
    load('../data/full_sampling_loc.mat');
    
    output_file = '../result/gamp_rec_3d';     % the prefix for the 3d output 

    % wavelet transform paramters
    nlevel = 4;     % wavelet transform level, usually 3-4 is enough
    wave_idx = 6;   % specify db1 to db8, use db6 to balance complexity and performance

    % least square via conjungate gradient paramters
    cg_tol=1e-4;            % cg convergence threshold
    max_cg_ite = 20;        % the number of cg iterations to compute the least square solution of multi-echo images
    est_init_marker = 1;    % if 1, compute LSQ solution; if 0, the LSQ solution was already computed, so just read LSQ solution

    % gamp parameters
    cvg_thd = 1e-6;     % convergence threshold
    kappa = 1;          % damping rate for parameter estimation
    x_mag_max = 0.1;    % maximum value of magnitude image
    x_mag_min = eps;    % minimum value of magnitude image
    x0_max = 0.1;       % maximum value of the initial magnetization
    x0_min = eps;       % minimum value of the initial magnetization
    r2star_max = 0.25;  % maixmum value of r2star value
    r2star_min = 0;     % minimum value of r2star value

    tau_w_1 = 1e-12;    % the initial noise variance for the multi-echo image prior
    tau_w_2 = 1e-12;    % the initial noise variance for the r2star estimation
    
    damp_rate = [];     % the damping rate of AMP iterations (like the step size), it is a number between 0 and 1, 
    if (sampling_rate<=0.1)
        damp_rate = 0.5;    % when the sampling rate is low, use a smaller damping rate
    else 
        damp_rate = 1;      % when the sampling rate is high, use a larger dampling rate, "1" means no dampling at all
    end

    % AMP is more computationally efficient, it converges faster than L1 and thus requires less iterations
    max_pe_spar_ite = 50;  % the number of iterations used by AMP to compute multi-echo image prior
    max_pe_ite = 50;       % the number of iterations used by AMP to recover initial magnetization, T2star/R2star map and phase images
    max_pe_est_ite = 5;    % the number of inner iterations to estimate paramters

    % create dictionary for estimation of r2star
    r2_dict_val = linspace(0,r2star_max,1001);
    exp_dict = zeros(Ne, length(r2_dict_val));
    for (i=1:length(r2_dict_val))
        exp_dict(:,i) = exp(-echo_time*r2_dict_val(i));
    end
    exp_dict_sq_norm = sum(exp_dict.^2,1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 2: generate Poisson-disk pattern %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    loc_index = reshape(1:(sy*sz), [sy sz]);    % sampling location indices
    
    sampling_array = ceil(poissonDisc_complement([sy, sz], sampling_pd_interval, [], num_samples-ncalib*ncalib));
    sampling_array = unique(sampling_array, 'row');

    % make sure the AC region is fully sampled
    % assume both sy and sz are even numbers
    % and that ncalib is also an even number

    ac_array = zeros(ncalib*ncalib, 2);
    ac_idx = 1;
    for (ac_y = (sy/2-ncalib/2+1):(sy/2+ncalib/2))
        for (ac_z = (sz/2-ncalib/2+1):(sz/2+ncalib/2))
            ac_array(ac_idx, :) = [ac_y ac_z];
            ac_idx = ac_idx+1;
        end
    end

    sampling_array = [sampling_array; ac_array];
    sampling_array = unique(sampling_array, 'row');

    sampling_loc = [];
    for (j=1:size(sampling_array,1))
        if (full_sampling_loc(sampling_array(j,1), sampling_array(j,2))==1)
            sampling_loc = [sampling_loc loc_index(sampling_array(j,1), sampling_array(j,2))];    % some times the sampling_loc is less than num_samples
        end
    end
    sampling_loc = unique(sampling_loc);    % only keep unique locations

    % if the number of samples is still less than it should be, poissonDisk could not help. Use random sampling then
    if (length(sampling_loc)<num_samples)
        left_locations = 1:(sy*sz);
        left_locations(full_sampling_loc(:)==0)=0;
        left_locations(sampling_loc)=0; 
        left_locations(left_locations==0)=[];
        sampling_loc = [sampling_loc randsample(left_locations, num_samples-length(sampling_loc))];
    end
    
    sampling_vect = zeros(num_samples, Ne);     % sampling vector
    for (echo_idx = 1:Ne)
        sampling_vect(:,echo_idx) = sort(sampling_loc);
    end 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 2.1: read data into noisy_measure %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    noisy_measure = zeros(sx, num_samples, Ne, Nc);
    for (echo_idx = 1:Ne)

        load(strcat('../data/echo_', num2str(echo_idx), '_3d.mat'))

        % generate measurements
        for (j=1:Nc)
            for (k=1:sx)
                X_fft_noisy_tmp = squeeze(data(k,:,:,j));
                noisy_measure(k,:,echo_idx,j) = X_fft_noisy_tmp(sampling_vect(:,echo_idx));
            end
        end

        % clear data
        clear data;
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% read the sensitivity 3D maps estiamted via espirit %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('../data/sensitivity_map_3d.mat')  % get maps_3d


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% use wavelet transform to enforce the sparsity of wavelet coefficients %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X0 = zeros(sx, sy, sz);

    % construct the sensing matrix
    dwtmode('per');
    C1=wavedec3(X0,nlevel,'db1'); 
    ncoef1=length(C1.dec);
    C2=wavedec3(X0,nlevel,'db2'); 
    ncoef2=length(C2.dec);
    C3=wavedec3(X0,nlevel,'db3'); 
    ncoef3=length(C3.dec);
    C4=wavedec3(X0,nlevel,'db4'); 
    ncoef4=length(C4.dec);
    C5=wavedec3(X0,nlevel,'db5'); 
    ncoef5=length(C5.dec);
    C6=wavedec3(X0,nlevel,'db6'); 
    ncoef6=length(C6.dec);
    C7=wavedec3(X0,nlevel,'db7'); 
    ncoef7=length(C7.dec);
    C8=wavedec3(X0,nlevel,'db8'); 
    ncoef8=length(C8.dec);


    switch wave_idx
    case 1
        Psi = @(x) [wavedec3(x,nlevel,'db1')']; 
        Psit = @(x) (waverec3(x));
        wav_vessel = C1;
    case 2
        Psi = @(x) [wavedec3(x,nlevel,'db2')']; 
        Psit = @(x) (waverec3(x));
        wav_vessel = C2;
    case 3
        Psi = @(x) [wavedec3(x,nlevel,'db3')']; 
        Psit = @(x) (waverec3(x));
        wav_vessel = C3;
    case 4
        Psi = @(x) [wavedec3(x,nlevel,'db4')']; 
        Psit = @(x) (waverec3(x));
        wav_vessel = C4;
    case 5
        Psi = @(x) [wavedec3(x,nlevel,'db5')']; 
        Psit = @(x) (waverec3(x));
        wav_vessel = C5;
    case 6
        Psi = @(x) [wavedec3(x,nlevel,'db6')']; 
        Psit = @(x) (waverec3(x));
        wav_vessel = C6;
    case 7
        Psi = @(x) [wavedec3(x,nlevel,'db7')']; 
        Psit = @(x) (waverec3(x));
        wav_vessel = C7;
    otherwise
        Psi = @(x) [wavedec3(x,nlevel,'db8')']; 
        Psit = @(x) (waverec3(x));
        wav_vessel = C8;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 3: use conjungate gradient descent to find %%
    %% the least squares solution with minimum l2 norm %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    y_echo = noisy_measure;
    clear noisy_measure;

    mat_sz = [sx sy sz];

    if (est_init_marker) 
        X_init = zeros(sx,sy,sz,Ne);
        for (j=1:Ne)
            d_cg = zeros(sx, sy, sz);
            X_init_tmp = zeros(sx, sy, sz);
            for (k=1:Nc)
                d_cg = d_cg - sx*sy*sz*conj(maps_3d(:,:,:,k)).*At_op_3d_cylinder(A_op_3d_cylinder(maps_3d(:,:,:,k).*X_init_tmp, sampling_vect(:,j))-y_echo(:,:,j,k), sampling_vect(:,j), mat_sz);
            end
            r_cg = d_cg;
            for (ite=1:max_cg_ite)  
                a_cg_n = sum(conj(r_cg).*r_cg, 'all');
                a_cg_d = 0;
                for (k=1:Nc)
                    a_cg_d = a_cg_d + sum(sx*sy*sz*(conj(d_cg).*conj(maps_3d(:,:,:,k)).*At_op_3d_cylinder(A_op_3d_cylinder(maps_3d(:,:,:,k).*d_cg, sampling_vect(:,j)), sampling_vect(:,j), mat_sz)), 'all');
                end
                a_cg_d = real(a_cg_d);
                a_cg = a_cg_n / a_cg_d;
                X_init_tmp_pre = X_init_tmp;
                X_init_tmp = X_init_tmp + a_cg * d_cg;
                cvg_cg_val = norm(X_init_tmp(:)-X_init_tmp_pre(:), 'fro')/norm(X_init_tmp(:), 'fro');
                fprintf('Echo %d, cvg val %d\n', j, cvg_cg_val)
                if (cvg_cg_val<cg_tol)
                    break;
                end
                r_cg_new = r_cg;
                for (k=1:Nc)
                    r_cg_new = r_cg_new - a_cg*sx*sy*sz*conj(maps_3d(:,:,:,k)).*At_op_3d_cylinder(A_op_3d_cylinder(maps_3d(:,:,:,k).*d_cg, sampling_vect(:,j)), sampling_vect(:,j), mat_sz);
                end
                b_cg = sum(conj(r_cg_new).*r_cg_new, 'all')/sum(conj(r_cg).*r_cg, 'all');
                d_cg = r_cg_new + b_cg*d_cg;
                r_cg = r_cg_new;
            end
            X_init(:,:,:,j) = X_init_tmp;
        end 

        save(strcat(output_file, '_X_init'), 'X_init')

    else

        load(strcat(output_file, '_X_init'))

    end 


    % compute initializations for the initial magnetization H_exp and R2* map 
    Et = abs(X_init);
    Et(Et<eps) = eps;
    Et_log = log(Et);

    warning('off','all')
    H = zeros(sx,sy,sz);
    R2_star = zeros(sx,sy,sz);
    for (idx_x=1:sx)
        fprintf('%d\n',idx_x)
        for (idx_y=1:sy)
            for (idx_z=1:sz)
            T_mat = [repmat(1, Ne, 1) -echo_time];
            Et_tmp = squeeze(Et(idx_x,idx_y,idx_z,:));
            T_mat(:,1) = T_mat(:,1).*Et_tmp;
            T_mat(:,2) = T_mat(:,2).*Et_tmp;
            Et_log_tmp = squeeze(Et_log(idx_x,idx_y,idx_z,:));
            Et_log_tmp = Et_log_tmp.*Et_tmp;
            T_mat=inv(T_mat'*T_mat)*T_mat';
            H_R2_star_seq = T_mat*Et_log_tmp;
            H(idx_x,idx_y,idx_z)=H_R2_star_seq(1);
            R2_star(idx_x,idx_y,idx_z)=H_R2_star_seq(2);
            end
        end
    end
    warning('on','all')

    H_exp = exp(H);
    H_exp = min(H_exp, x0_max);
    H_exp = max(H_exp, x0_min);
    
    H_exp_psi_struct = Psi(H_exp);
    H_exp_psi = extract_3d_wav_coef(H_exp_psi_struct);
    
    R2_star = min(R2_star, r2star_max);
    R2_star = max(R2_star, r2star_min);

    wav_coef_len = length(H_exp_psi);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 4: gamp reconstruction %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	% initialize the gamp parameters
    gamp_par.damp_rate = damp_rate;
    gamp_par.max_pe_ite = max_pe_ite;
    gamp_par.max_pe_spar_ite = max_pe_spar_ite;
    gamp_par.max_pe_est_ite = max_pe_est_ite;

    gamp_par.cvg_thd = cvg_thd;
    gamp_par.kappa = kappa;
    gamp_par.echo_time = echo_time;
    gamp_par.sx = sx;
    gamp_par.sy = sy;
    gamp_par.sz = sz;
    gamp_par.Ne = Ne;
    gamp_par.Nc = Nc;

    gamp_par.x_mag_max = x_mag_max;
    gamp_par.x_mag_min = x_mag_min;
    gamp_par.x0_max = x0_max;
    gamp_par.x0_min = x0_min;
    gamp_par.r2star_max = r2star_max;
    gamp_par.r2star_min = r2star_min;

	
    % initialize the variables
    x_hat_meas = X_init;
    %clear X_init;

    % initialize gamp variables
    x_hat_meas_psi = zeros([wav_coef_len Ne]);
    for (i=1:Ne)
        x_hat_meas_psi_struct_tmp = Psi(x_hat_meas(:,:,:,i));
        x_hat_meas_psi(:,i) = extract_3d_wav_coef(x_hat_meas_psi_struct_tmp);
    end
    tau_x_meas_psi = var(x_hat_meas_psi(:));

    gamp_par.tau_x_meas_psi = tau_x_meas_psi;
    gamp_par.x_hat_meas_psi = x_hat_meas_psi;
    gamp_par.s_hat_meas_1 = zeros(size(y_echo));
    gamp_par.s_hat_meas_2 = zeros(size(y_echo));

    % initialize the distribution parameters
    lambda_x_hat_psi = [];
    for (i = 1:Ne)
        x_hat_psi_tmp = x_hat_meas_psi(:,i);
        lambda_x_hat_psi = [lambda_x_hat_psi; 1/sqrt(var(abs(x_hat_psi_tmp))/2)];
    end
    input_par.lambda_x_hat_psi = lambda_x_hat_psi;

    input_par.lambda_x0_psi = 1/sqrt(var(H_exp_psi)/2);

    % output variance parameter
    output_par.tau_w_1 = tau_w_1;
    output_par.tau_w_2 = tau_w_2;
	
    % construct measurement operator needed for multi-echo reconstruction 
    hS = @(in1, in2) A_op_3d_cylinder(in1, in2);
    hSt = @(in1, in2) At_op_3d_cylinder(in1, in2, mat_sz);

    E_coef = @(in) extract_3d_wav_coef(in);
    C_struct = @(in) construct_3d_wav_struct(in, wav_vessel);

    M=num_samples*sx*Ne*Nc;
    N=wav_coef_len*Ne;
    A_2_echo_3d = A_2_echo_3d_cylinder_LinTrans(M,N,maps_3d,mat_sz,num_samples,Ne,Nc,wav_coef_len,sampling_vect,Psit, Psi,C_struct,E_coef,hS,hSt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% multi-echo reconstruction %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [res, input_par_new, output_par_new] = gamp_mri_multi_echo_3d_nw(A_2_echo_3d, y_echo, gamp_par, input_par, output_par);

    save(strcat(output_file, '_multi_echo_res'), 'res', '-v7.3')
    save(strcat(output_file, '_multi_echo_input_par_new'), 'input_par_new')
    save(strcat(output_file, '_multi_echo_output_par_new'), 'output_par_new')

    load(strcat(output_file, '_multi_echo_res'))
    load(strcat(output_file, '_multi_echo_input_par_new'))
    load(strcat(output_file, '_multi_echo_output_par_new'))


    % update the parameters and variables computed from multi-echo reconstruction
    gamp_par.x_hat_meas_psi = res.x_hat_meas_psi;
    gamp_par.tau_x_meas_psi = res.tau_x_meas_psi;
    gamp_par.s_hat_meas_1 = res.s_hat_meas_1;

    input_par.lambda_x_hat_psi = input_par_new.lambda_x_hat_psi;
    output_par.tau_w_1 = output_par_new.tau_w_1;

    %clear res;
    %clear A_2_echo_3d;

    % construct measurement operators needed for r2star reconstruciton
    M=sx*sy*sz*Ne;
    N=wav_coef_len;
    A_x0_r2star_3d = A_x0_r2star_3d_LinTrans(M,N,mat_sz,wav_coef_len,Ne,echo_time,Psit,Psi,C_struct,E_coef);

    M=sx*num_samples*Ne*Nc;
    N=sx*sy*sz*Ne;
    A_2_echo_nw_3d = A_2_echo_nw_3d_cylinder_LinTrans(M,N,maps_3d,mat_sz,num_samples,Ne,Nc,sampling_vect,hS,hSt);

    M=sx*sy*sz*Ne;
    N=wav_coef_len*Ne;
    A_wav_3d = A_wav_3d_LinTrans(M,N,mat_sz,wav_coef_len,Ne,Psit,Psi,C_struct,E_coef);

    M=sx*sy*sz;
    N=wav_coef_len;
    A_wav_single_3d = A_wav_single_3d_LinTrans(M,N,mat_sz,wav_coef_len,Psit,Psi,C_struct,E_coef);

    % initialize the exp part from the least square solution
    gamp_par.x_hat_init = x_hat_meas;

    % initialize dictionary needed for estimating r2star values
    gamp_par.r2_dict_val = r2_dict_val;
    gamp_par.exp_dict = exp_dict;
    gamp_par.exp_dict_sq_norm = exp_dict_sq_norm;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% r2star reconstruciton %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    [res, input_par_new, output_par_new] = gamp_mri_x0_r2star_3d_nw_dict(A_2_echo_nw_3d, A_wav_3d, A_wav_single_3d, A_x0_r2star_3d, y_echo, gamp_par, input_par, output_par);

    save(strcat(output_file, '_res'), 'res', '-v7.3')
    save(strcat(output_file, '_input_par_new'), 'input_par_new')
    save(strcat(output_file, '_output_par_new'), 'output_par_new')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot recovered images %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % note that depending on the dataset, ifftshift or fftshift might be needed
    x0_rec = ifftshift(Psit(C_struct(res.x0_psi)));    % recovered initial magnetization
    r2star_rec = ifftshift(res.r2star);                % recovered R2star map
    x_hat_rec = res.x_hat_all;              % recovered multi-echo images
    for (i=1:Ne)
        x_hat_rec(:,:,:,i) = ifftshift(x_hat_rec(:,:,:,i));
    end

    % note that the display range needs to be set properly
    figure; imshow3D(x0_rec,[0,0.015])
    figure; imshow3D(r2star_rec,[0,0.05])
    figure; imshow3D(real(x_hat_rec(:,:,:,1)),[-0.015,0.015])
