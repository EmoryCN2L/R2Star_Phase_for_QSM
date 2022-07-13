    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 1: Setting up parameters %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Regularization parameters lambda_1 and lambda_2 need to be *tuned* for different types of datasets
    lambda_1 = 100;     % regularization parameter for recovering multi-echo images prior
    lambda_2 = 5e-5;    % regularization parameter for recovering initial magnetization, T2* maps

    % ESPIRiT parameters
    ncalib = 24;        % size of central k-space, use 24 calibration lines to estimate senstivity map, make sure to do full sampling in the AC region

    % dataset parameters
    sy = 204;       % size along y-direction
    sz = 164;       % size along z-direction
    Nc = 16;        % the number of channels (coils)
    Ne = 4;         % the number of echoes
    echo_time = [4 12 20 28].'; % echo time in ms

    LN_num = 10;                        % computing thread
    LN = maxNumCompThreads( LN_num );   % set the largest number of computing threads

    % generate Poisson sampling pattern accroding to poisson disk
    sampling_rate = 0.15;        % make sure to do full sampling in the AC region
    sampling_pd_interval = 2;   % the minimum distance between any two locations, needs   to be tuned for best performance
    num_samples = round(sy*sz*sampling_rate);   % the number of samples

    % for this dataset, sx direction is fully sampled, undersampling takes place in sy-sz plane
    % full_sampling_loc should be a mask of size sy by sz, the mask contains 0-1 values 
    load('../data/Sim1/full_sampling_loc_204_164.mat');
    load('../data/Sim1/mask_2d.mat')    % mask for the brain region

    output_file = '../result/Sim1_l1_reg_rec_2d';     % the prefix for the 3d output 

    % wavelet transform paramters
    nlevel = 4;     % wavelet transform level, usually 3-4 is enough
    wave_idx = 6;   % specify db1 to db8, use db6 to balance complexity and performance

    % least square via conjungate gradient paramters
    cg_tol=1e-4;            % cg convergence threshold
    max_cg_ite = 20;        % the number of cg iterations to compute the least square solution of multi-echo images
    est_init_marker = 1;    % if 1, compute LSQ solution; if 0, the LSQ solution was      already computed, so just read LSQ solution

    % recovery parameters
    cvg_thd = 1e-6;     % convergence threshold
    x_mag_max = 500;    % maximum value of magnitude image ** needs to be updated accordingly **
    x_mag_min = eps;    % minimum value of magnitude image
    x0_max = 500;       % maximum value of the initial magnetization ** needs to be updated accordingly **
    x0_min = eps;       % minimum value of the initial magnetization
    r2star_max = 0.25;  % maixmum value of r2star value
    r2star_min = 0;     % minimum value of r2star value

    maxiter = 100;      % the number of FISTA iterations
    maxinneriter = 10;  % the number of inner iterations to recover the initial magnetization


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

    loc_index = reshape(1:(sy*sz), [sy sz]);

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
            sampling_loc = [sampling_loc loc_index(sampling_array(j,1), sampling_array(j, 2))];    % some times the sampling_loc is less than num_samples
        end
    end
    sampling_loc = unique(sampling_loc);    % only keep unique locations

    % if the number of samples is still less than it should be, poissonDisk could not     help. Use random sampling then
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

    noisy_measure = zeros(num_samples, Ne, Nc);

    load(strcat('../data/Sim1/slice_2d.mat'))
    for (i=1:Ne)
        for (j=1:Nc)
            data_tmp = squeeze(data(:,:,i,j));
            noisy_measure(:,i,j) = data_tmp(sampling_vect(:,i));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the sensitivity 3D maps estiamted via espirit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('../data/Sim1/sensitivity_map_2d.mat')  % get maps


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use wavelet transform to enforce the sparsity of wavelet coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X0 = zeros(sy, sz);

    % construct the sensing matrix
    dwtmode('per');
    [C1,S1]=wavedec2(X0,nlevel,'db1');
    ncoef1=length(C1);
    [C2,S2]=wavedec2(X0,nlevel,'db2');
    ncoef2=length(C2);
    [C3,S3]=wavedec2(X0,nlevel,'db3');
    ncoef3=length(C3);
    [C4,S4]=wavedec2(X0,nlevel,'db4');
    ncoef4=length(C4);
    [C5,S5]=wavedec2(X0,nlevel,'db5');
    ncoef5=length(C5);
    [C6,S6]=wavedec2(X0,nlevel,'db6');
    ncoef6=length(C6);
    [C7,S7]=wavedec2(X0,nlevel,'db7');
    ncoef7=length(C7);
    [C8,S8]=wavedec2(X0,nlevel,'db8');
    ncoef8=length(C8);


    switch wave_idx
    case 1
        Psi = @(x) [wavedec2(x,nlevel,'db1')']; 
        Psit = @(x) (waverec2(x(1:ncoef1),S1,'db1'));
    case 2
        Psi = @(x) [wavedec2(x,nlevel,'db2')']; 
        Psit = @(x) (waverec2(x(1:ncoef2),S2,'db2'));
    case 3
        Psi = @(x) [wavedec2(x,nlevel,'db3')']; 
        Psit = @(x) (waverec2(x(1:ncoef3),S3,'db3'));
    case 4
        Psi = @(x) [wavedec2(x,nlevel,'db4')']; 
        Psit = @(x) (waverec2(x(1:ncoef4),S4,'db4'));
    case 5
        Psi = @(x) [wavedec2(x,nlevel,'db5')']; 
        Psit = @(x) (waverec2(x(1:ncoef5),S5,'db5'));
    case 6
        Psi = @(x) [wavedec2(x,nlevel,'db6')']; 
        Psit = @(x) (waverec2(x(1:ncoef6),S6,'db6'));
    case 7
        Psi = @(x) [wavedec2(x,nlevel,'db7')']; 
        Psit = @(x) (waverec2(x(1:ncoef7),S7,'db7'));
    otherwise
        Psi = @(x) [wavedec2(x,nlevel,'db8')']; 
        Psit = @(x) (waverec2(x(1:ncoef8),S8,'db8'));
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 3: use conjungate gradient descent to find %%
    %% the least squares solution with minimum l2 norm %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y_echo = noisy_measure;
    clear noisy_measure;

    mat_sz = [sy sz];

    if (est_init_marker)
        X_init = zeros(sy,sz,Ne);
        for (j=1:Ne)
            d_cg = zeros(sy, sz);
            X_init_tmp = zeros(sy, sz);
            for (k=1:Nc)
                d_cg = d_cg - sy*sz*conj(maps(:,:,k)).*At_op(A_op(maps(:,:,k).*X_init_tmp, sampling_vect(:,j))-y_echo(:,j,k), sampling_vect(:,j), mat_sz);
            end
            r_cg = d_cg;
            for (ite=1:max_cg_ite)   % max iteration is 100
                a_cg_n = sum(conj(r_cg).*r_cg, 'all');
                a_cg_d = 0;
                for (k=1:Nc)
                    a_cg_d = a_cg_d + real(sum(sy*sz*(conj(d_cg).*conj(maps(:,:,k)).*At_op(A_op(maps(:,:,k).*d_cg, sampling_vect(:,j)), sampling_vect(:,j), mat_sz)), 'all'));
                end
                a_cg = a_cg_n / a_cg_d;
                X_init_tmp_pre = X_init_tmp;
                X_init_tmp = X_init_tmp + a_cg * d_cg;
                cvg_cg_val = norm(X_init_tmp-X_init_tmp_pre, 'fro')/norm(X_init_tmp, 'fro');
                fprintf('Echo %d, cvg val %d\n', j, cvg_cg_val)
                if (cvg_cg_val<cg_tol)
                    break;
                end
                r_cg_new = r_cg;
                for (k=1:Nc)
                    r_cg_new = r_cg_new - a_cg*sy*sz*conj(maps(:,:,k)).*At_op(A_op(maps(:,:,k).*d_cg, sampling_vect(:,j)), sampling_vect(:,j), mat_sz);
                end
                b_cg = sum(conj(r_cg_new).*r_cg_new, 'all')/sum(conj(r_cg).*r_cg, 'all');
                d_cg = r_cg_new + b_cg*d_cg;
                r_cg = r_cg_new;
            end
            X_init(:,:,j) = X_init_tmp;
        end

        save(strcat(output_file, '_X_init'), 'X_init')

    else

        load(strcat(output_file, '_X_init'))

    end


    par.X_init = X_init;    % initialize the estimated image with all zeros
    par.maxiter = maxiter;
    par.tol = cvg_thd;            % adjust this accordingly
    kappa = zeros(1, Nc);
    for (j=1:Nc)
        kappa(j) = 2*sy*sz*(max(abs(maps(:,:,j)), [], 'all')).^2;
    end 
    kappa = max(kappa(:));
    par.kappa = kappa;             % should be twice the largest eigenvalue of (At*A),     need to tune to make sure
    par.s_vect = sampling_vect;
    par.x_mag_max = x_mag_max;
    par.x_mag_min = x_mag_min;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 4: L1-min reconstruction %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [X_rec_echo, fun_val_block_1_cur_0, fun_val_block_1_cur_1] = recovery_multi_wave_l1_reg_2d(y_echo, maps, lambda_1, Psi, Psit, par);   

    save(strcat(output_file, '_multi_echo'), 'X_rec_echo')
    save(strcat(output_file, '_multi_echo_fun_val_0'), 'fun_val_block_1_cur_0')
    save(strcat(output_file, '_multi_echo_fun_val_1'), 'fun_val_block_1_cur_1')


    % compute initializations for the initial magnetization H_exp and R2* map
    Et = abs(X_rec_echo);
    Et(Et<eps) = eps;
    Et_log = log(Et);

    warning('off','all')
    H = zeros(sy,sz);
    R2_star = zeros(sy,sz);
    for (idx_y=1:sy)
        for (idx_z=1:sz)
            T_mat = [repmat(1, Ne, 1) -echo_time];
            Et_tmp = squeeze(Et(idx_y,idx_z,:));
            T_mat(:,1) = T_mat(:,1).*Et_tmp;
            T_mat(:,2) = T_mat(:,2).*Et_tmp;
            Et_log_tmp = squeeze(Et_log(idx_y,idx_z,:));
            Et_log_tmp = Et_log_tmp.*Et_tmp;
            T_mat=inv(T_mat'*T_mat)*T_mat';
            H_R2_star_seq = T_mat*Et_log_tmp;
            H(idx_y,idx_z)=H_R2_star_seq(1);
            R2_star(idx_y,idx_z)=H_R2_star_seq(2);
        end
    end 
    warning('on','all')

    H_exp = exp(H);
    H_exp = min(H_exp, x0_max);
    H_exp = max(H_exp, x0_min);

    H_exp_psi = Psi(H_exp);

    R2_star = min(R2_star, r2star_max);
    R2_star = max(R2_star, r2star_min);

    wav_coef_len = length(H_exp_psi);


    par.H_exp = H_exp;
    par.H_exp_max = x0_max;
    par.H_exp_min = x0_min; 
    par.V = R2_star;
    par.V_max = r2star_max;
    par.V_min = r2star_min;
    par.maxiter = maxiter;
    par.maxinneriter = maxinneriter;  % since the estimation of H_exp is continuous, we can set its iteration number smaller
    par.tol = cvg_thd;            % adjust this accordingly
    par.echo_time = echo_time;       % adjust this accordingly

    par.r2_dict_val = r2_dict_val;
    par.exp_dict = exp_dict;
    par.exp_dict_sq_norm = exp_dict_sq_norm;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reconstruct r2star and the initial magnetization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [X0_rec, R2_star_rec, fun_val_block_1_cur_0, fun_val_block_1_cur_1] = exponential_fit_2d_x0_inner_r2_l1_reg_dict(Et, lambda_2, Psi, Psit, par);

    save(strcat(output_file, '_X0'), 'X0_rec')
    save(strcat(output_file, '_R2_star'), 'R2_star_rec')
    save(strcat(output_file, '_x0_r2star_fun_val_0'), 'fun_val_block_1_cur_0')
    save(strcat(output_file, '_x0_r2star_fun_val_1'), 'fun_val_block_1_cur_1')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot recovered images %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % note that depending on the dataset, ifftshift or fftshift might be needed
    x0_rec = X0_rec;             % recovered initial magnetization
    r2star_rec = R2_star_rec;    % recovered R2star map
    x_hat_rec = X_rec_echo;                 % recovered multi-echo images

    % extract the brain with a brain mask
    x0_rec(mask_2d==0) = 0;
    r2star_rec(mask_2d==0) = 0;
    x_hat_rec(mask_2d==0) = 0;

    % note that the display range needs to be set properly
    figure; imshow(x0_rec,[0,200])
    figure; imshow(r2star_rec,[0,0.1])
    figure; imshow(real(x_hat_rec(:,:,1)),[-200,200])
