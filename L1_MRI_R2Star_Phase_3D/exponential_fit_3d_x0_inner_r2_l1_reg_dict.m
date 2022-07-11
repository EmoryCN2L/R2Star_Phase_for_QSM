function [H_exp, V, fun_val_block_1_cur_0, fun_val_block_1_cur_1]=exponential_fit_3d_x0_inner_r2_l1_reg_dict(Et, lambda_2, Psi, Psit, E_coef, C_struct, par)

    % regularization on X0
    
    % Et        : the recovered magnitude echo images of size sx by sy by Ne
    %			: sx is the resolution in the x direction
    %			: sy is the resolution in the y direction
    %			: Ne is the number of echoes
    
    % par       : various parameters
    
    % Assigning parameres according to par

    maxiter = par.maxiter;          % the maximum number of fista iterations
    maxinneriter = par.maxinneriter;    % the maximum number of inner iteration w.r.t. H_exp using fista
    tol     = par.tol;              % convergence criterion, can be set to 1e-6
    echo_time = par.echo_time;		% the echo times, should be scaled properly to speed up convergence
    

    % maximum and minimum value of Rt_abs

    sx = size(Et,1);
    sy = size(Et,2);
    sz = size(Et,3);
    echo_num = size(Et,4);

    H_exp_min = par.H_exp_min;
    H_exp_max = par.H_exp_max;
    V_min = par.V_min;  % or eps???
    V_max = par.V_max;

    
    % initialize H, H_exp, and V
    %H_exp = repmat(eps, [sx sy]);
    H_exp = par.H_exp;
    H_exp = max(H_exp, H_exp_min);

    %V = repmat(eps, [sx sy]);
    V = par.V;
    V = max(V, V_min);

    exp_dict = par.exp_dict;
    exp_dict_sq_norm = par.exp_dict_sq_norm;
    r2_dict_val = par.r2_dict_val;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % admm block 2: optimize with respect to H and V %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k_t = 1;
    k_tm1 = 1;
    St_H_exp = H_exp;
    St_V = V;
    Stm1_H_exp = St_H_exp;

    fun_val_block_1_cur = 0;
    cvg_count_block_1 = 0;
    cvg_count_block_2 = 0;
    for (ii=1:maxiter)

        for (jj=1:maxinneriter)
        H_exp = St_H_exp + ((k_tm1-1)/k_t)*(St_H_exp-Stm1_H_exp);

        H_exp = min(H_exp, H_exp_max);
        H_exp = max(H_exp, H_exp_min);

        % optimize with respect to H_exp
        Lpz_H_exp = zeros(echo_num,1);
        T = zeros(size(Et));
        for (j=1:echo_num)
            T(:,:,:,j) = exp(-echo_time(j)*V);
            Lpz_H_exp(j) = 2*max(T(:,:,:,j).^2,[],'all');
        end

        Lpz_H_exp = repmat(max(Lpz_H_exp), size(Lpz_H_exp));

        H_exp_par_n = zeros(sx,sy,sz);
        H_exp_par_d = 0;
        for (j=1:echo_num)
            H_exp_par_n = H_exp_par_n + Lpz_H_exp(j)*H_exp-2*(H_exp.*T(:,:,:,j)-Et(:,:,:,j)).*T(:,:,:,j);
            H_exp_par_d = H_exp_par_d + Lpz_H_exp(j);
        end

        H_exp = H_exp_par_n/H_exp_par_d;
        H_exp_Psi = Psi(H_exp);
        H_exp_Psi = extract_3d_wav_coef(H_exp_Psi);
        H_exp_Psi = max(abs(H_exp_Psi)-lambda_2/H_exp_par_d, 0) .* sign(H_exp_Psi);
        H_exp_Psi = C_struct(H_exp_Psi);
        H_exp = Psit(H_exp_Psi);
        H_exp = min(H_exp, H_exp_max);
        H_exp = max(H_exp, H_exp_min);

        Stp1_H_exp = H_exp;

        fun_val_block_1_pre = fun_val_block_1_cur;
        fun_val_block_1_cur_0 = 0;
        for (j=1:echo_num)
            fun_val_block_1_cur_0 = fun_val_block_1_cur_0 + sum((Et(:,:,:,j)-H_exp.*exp(-echo_time(j)*V)).^2, 'all');
        end
        fun_val_block_1_cur_1 = sum(abs(extract_3d_wav_coef(Psi(H_exp))));
        fun_val_block_1_cur = fun_val_block_1_cur_0 + lambda_2 * fun_val_block_1_cur_1;
        %cvg_val_block_1 = abs((fun_val_block_1_pre-fun_val_block_1_cur)/fun_val_block_1_cur);

        cvg_val_block_H_exp_1 = sqrt(sum(abs(Stp1_H_exp - St_H_exp).^2, 'all'))/sqrt(sum(abs(Stp1_H_exp).^2, 'all'));
        cvg_val_block_1 = max(cvg_val_block_H_exp_1);

        if (cvg_val_block_1<tol)
            cvg_count_block_1 = cvg_count_block_1 + 1;
        end

        if (cvg_count_block_1>10)
            break;
        end

        fprintf('Block 2: %d   %d   %d\n', jj, fun_val_block_1_cur, cvg_val_block_H_exp_1)

        % update
        k_tp1 = 0.5*(1+sqrt(1+4*k_t*k_t)) ;

        k_tm1 = k_t ;
        k_t = k_tp1 ;

        Stm1_H_exp = St_H_exp;
        St_H_exp = Stp1_H_exp;
        end

        % optimize with respect to R2

        V = exp_lsq_fit_3d_r2_componentwise_dict(Et, H_exp, exp_dict, exp_dict_sq_norm, r2_dict_val);
        V = min(V, V_max);
        V = max(V, V_min);

        Stp1_V = V;
        cvg_val_block_V_2 = sqrt(sum(abs(Stp1_V - St_V).^2, 'all'))/sqrt(sum(abs(Stp1_V).^2, 'all'));
        fprintf('V opt: %d  %d\n', ii, cvg_val_block_V_2);

	if (cvg_val_block_V_2<tol)
		cvg_count_block_2 = cvg_count_block_2 + 1;
	end

	if (cvg_count_block_2>10)
		break;
	end

        St_V = Stp1_V;

    end

end
