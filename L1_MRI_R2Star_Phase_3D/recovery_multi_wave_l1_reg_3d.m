function [Stp1, fun_val_block_1_cur_0, fun_val_block_1_cur_1] = recovery_multi_wave_l1_reg_3d(Y, maps, lambda_1, Psi, Psit, E_coef, C_struct, par)

    % adapted from previous code - S.H. 10/09/2019
    % optimize with respect to magnitude first and phase
    % use bisection to find the global optimum
    % add optimization of phase
    % with inner iterations
    % with real fista acceleration
    
    % Image recovery from linear measurements via sparsity averaging: 
    %
    % R. E. Carrillo, et al. “Sparsity averaging for compressive imaging,” IEEE Signal Processing Letters, vol. 20, no. 6, pp. 591–594, June 2013.
    % S. Setzer, Split Bregman algorithm, Douglas-Rachford splitting and frame shrinkage, pp. 464–476, Springer Berlin Heidelberg, Berlin, Heidelberg, June 2009.
    % 
    % By Shuai Huang, The Johns Hopkins University
    % Email: shuang40@jhu.edu 
    % Date: 12/20/2018

    % Y         : Observation/measurement
    % A         : structual random matrix operator
    % At        : transpose of A
    % lambda    : the regularization parameter
    % Psi       : overcomplete wavelet basis operator
    % Psit      : inverse overcomplete wavelet transform operator
    % Q         : a term used by Alternating Split Bregman Shrinkage (ASBS) algorithm
    % par       : various parameters
    
    % St        : the recovered images
    % X         : the sparse wavelet coefficient
    
    % Assigning parameres according to par

    maxiter = par.maxiter;          % the maximum number of fista iterations
    St      = par.X_init;               % initialization of estimated image
    kappa   = par.kappa;            % the Lipschitz constant
    tol     = par.tol;              % convergence criterion, can be set to 1e-6

    
    s_vect  = par.s_vect;       % sampling vector         

    % maximum and minimum value of Rt_abs
    max_val = 1;
    min_val = eps;

    mat_sz = [size(St, 1) size(St, 2) size(St, 3)];
    sx = size(St,1);
    sy = size(St,2);
    sz = size(St,3);
    echo_num = size(St,4);
    coil_num = size(maps, 4);
    
    Rt = St;
    P  = zeros([size(St) echo_num coil_num]);   % proximal regularization
    
    % is this the right way to initialize???

    fun_val_cur=0;
    fun_val_seq = [];
    
    con_val = [];
    cvg_count = 0;

    Rt_abs = abs(St);   % should we use this one

    St_abs = abs(St);
    St_abs(St_abs<eps) = eps;
    Rt_pr = St./St_abs;

    Rt = Rt_abs .* Rt_pr;

    Rt_abs_max = par.x_mag_max;
    Rt_abs_min = par.x_mag_min;

    maps_conj = conj(maps);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % admm block 0: use L1 to initialize H, V, E %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % use FISTA to solve it
    % should we do fista on Rt or Rt_abs? let us do it on Rt first
    k_t = 1;
    k_tm1 = 1;
    St = Rt;
    Stm1 = St;

    fun_val_block_1_cur = 0;
    cvg_count_block_1 = 0;

    fun_val_block_1_cur_0 = 0;
    fun_val_block_1_cur_1 = zeros(1,echo_num);
    for (ii=1:maxiter)
    Rt = St + ((k_tm1-1)/k_t)*(St-Stm1);

    % optimize the echo image reconstruction
    for (j=1:echo_num)
        for (k=1:coil_num)
            % multiply by 2 here, should we do this for complex case as well?
            P(:,:,:,j,k)=Rt(:,:,:,j)-1/kappa*2*sx*sy*sz*maps_conj(:,:,:,k).*At_op_3d_cylinder(A_op_3d_cylinder(maps(:,:,:,k).*Rt(:,:,:,j), s_vect(:,j))-Y(:,:,j,k), s_vect(:,j), mat_sz);
        end
    end

    % Step 1: optimize over the echo images

    for (j=1:echo_num)
        Rt_tmp = zeros(sx,sy,sz);
        for (k=1:coil_num)
            Rt_tmp = Rt_tmp + P(:,:,:,j,k);
        end
        Rt_tmp = Rt_tmp/coil_num;
        normalized_tmp = kappa*coil_num;

        Rt_tmp_psi = Psi(Rt_tmp);
        Rt_tmp_psi = extract_3d_wav_coef(Rt_tmp_psi);
        Rt_tmp_psi = max(abs(Rt_tmp_psi)-lambda_1/normalized_tmp, 0).*sign(Rt_tmp_psi);
        Rt_tmp_psi = C_struct(Rt_tmp_psi);
        Rt_tmp = Psit(Rt_tmp_psi);
        Rt_abs_tmp = abs(Rt_tmp);
        Rt_sign_tmp = sign(Rt_tmp);
        Rt_abs_tmp = min(Rt_abs_tmp, Rt_abs_max);
        Rt_abs_tmp = max(Rt_abs_tmp, Rt_abs_min);
        Rt_tmp = Rt_abs_tmp .* Rt_sign_tmp;
        Rt(:,:,:,j) = Rt_tmp;
    end
    Stp1 = Rt;

    fun_val_block_1_pre = fun_val_block_1_cur;
    fun_val_block_1_cur_0 = 0;
    for (j=1:echo_num)
        for (k=1:coil_num)
            fun_val_block_1_cur_0 = fun_val_block_1_cur_0 + norm(Y(:,:,j,k)-A_op_3d_cylinder(maps(:,:,:,k).*Stp1(:,:,:,j), s_vect(:,j)), 'fro')^2;
        end
    end

    fun_val_block_1_cur_1 = zeros(1,echo_num);
    for (j=1:echo_num)
        fun_val_block_1_cur_1(j) = sum(abs(extract_3d_wav_coef(Psi(Stp1(:,:,:,j)))), 'all');
    end

    fun_val_block_1_cur = fun_val_block_1_cur_0 + lambda_1*sum(fun_val_block_1_cur_1);
    %cvg_val_block_1 = abs((fun_val_block_1_pre-fun_val_block_1_cur)/fun_val_block_1_cur);
    cvg_val_block_1 = sqrt(sum(abs(Stp1 - St).^2, 'all'))/sqrt(sum(abs(Stp1).^2, 'all'));

    if (cvg_val_block_1<tol)
        cvg_count_block_1 = cvg_count_block_1+1;
    end
    if (cvg_count_block_1>10)
        break;
    end

    fprintf('Block 0: %3d   %d   %d\n', ii, fun_val_block_1_cur, cvg_val_block_1)

    % update 
    k_tp1 = 0.5*(1+sqrt(1+4*k_t*k_t)) ;

    k_tm1 = k_t ;
    k_t = k_tp1 ;

    Stm1 = St ;
    St = Stp1 ;

    end
    
end
