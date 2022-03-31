function [res, input_par, output_par] = gamp_mri_x0_r2star_2d_nw_dict(A_2_echo_nw, A_wav, A_wav_single, A_x0_r2star, y, gamp_par, input_par, output_par)

    % no wavelet in A_2_echo_nw
    % did not embed parameter estimation, move estimation of x_hat_all to the last step
    % now embed parameter estimation

    % seems that i should not subtract info from measurements when passing info from exp to meas, because it will diverge, so in this case i am not subtracting, use new way to compute x_hat_exp and tau_x_exp

    % still diverges, move x_hat_all ahead
    % use damping on both blocks

    % move x_hat_all behind, do not apply regularization on x0_psi, no damping

    % use a bi block message passing one end is x_hat_all 

    % only pass magnitude info to exp block, the variance is automatically estimated, hope this will stop the rescontruction from becoming low resolution

    % enforce positive constraint on the r2star
	% set GAMP parameters
	max_pe_ite = gamp_par.max_pe_ite;
    max_pe_spar_ite = gamp_par.max_pe_spar_ite;
    max_pe_est_ite = gamp_par.max_pe_est_ite;
	cvg_thd = gamp_par.cvg_thd;
	kappa = gamp_par.kappa;
    damp_rate = gamp_par.damp_rate;

    tau_x_meas_psi = gamp_par.tau_x_meas_psi;
    x_hat_meas_psi = gamp_par.x_hat_meas_psi;
    s_hat_meas_1 = gamp_par.s_hat_meas_1;
    s_hat_meas_2 = gamp_par.s_hat_meas_2;

    x_hat_init = gamp_par.x_hat_init;

    % set input distribution parameters
    lambda_x0_psi  = input_par.lambda_x0_psi; % Bernoulli parameter
    lambda_x_hat_psi = input_par.lambda_x_hat_psi;

    tau_w_1 = output_par.tau_w_1;
    tau_w_2 = output_par.tau_w_2;

    echo_time = gamp_par.echo_time;
    sy = gamp_par.sy;
    sz = gamp_par.sz;
    Ne = gamp_par.Ne;
    Nc = gamp_par.Nc;

    x_mag_max = gamp_par.x_mag_max;
    x_mag_min = gamp_par.x_mag_min;
    x0_max = gamp_par.x0_max;
    x0_min = gamp_par.x0_min;
    r2star_max = gamp_par.r2star_max;
    r2star_min = gamp_par.r2star_min;

    exp_dict = gamp_par.exp_dict;
    exp_dict_sq_norm = gamp_par.exp_dict_sq_norm;
    r2_dict_val = gamp_par.r2_dict_val;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% initialize tau_w_exp with tau_x_mag_meas from the warm up step %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% first finish the sparse reconstruction %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tau_x_mag_meas_1 = A_wav.multSq(tau_x_meas_psi);
    x_mag_meas_1_complex = A_wav.mult(x_hat_meas_psi);

    %% compute the estimation of X0 and R2star
    %Et = abs(x_mag_meas_1_complex);
    Et = abs(x_hat_init);
    Et(Et<x0_min) = x0_min;
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

    H_exp_psi = A_wav_single.multTr(H_exp);

    R2_star = min(R2_star, r2star_max);
    R2_star = max(R2_star, r2star_min); 

    x0_psi = H_exp_psi;
    tau_x0_psi = var(x0_psi(:)) ;
    r2star = R2_star;

    % initialize the A_x0_r2star operator
    set(A_x0_r2star, 'R2star', r2star);
    A_x0_r2star.estFrob();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% warm up the iterations %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tau_p_exp = A_x0_r2star.multSq(tau_x0_psi);
    p_hat_exp = A_x0_r2star.mult(x0_psi); % - tau_p_exp .* s_hat_exp;

    tau_x_mag_exp = tau_p_exp;
    x_mag_exp = p_hat_exp;
    x_mag_exp = max(x_mag_exp, x_mag_min);  % make sure x_mag_exp is positive

    tau_x_mag_meas = tau_x_mag_meas_1;
    x_mag_meas_complex = x_mag_meas_1_complex;
    x_mag_meas = abs(x_mag_meas_complex);
    x_mag_meas = max(x_mag_meas, x_mag_min);

    tau_x_mag_all = (tau_x_mag_exp*tau_x_mag_meas) / (tau_x_mag_exp+tau_x_mag_meas);
    x_mag_all = (tau_x_mag_exp*x_mag_meas+tau_x_mag_meas*x_mag_exp) / (tau_x_mag_exp+tau_x_mag_meas);

    x_hat_all = x_mag_all .* sign(x_mag_meas_complex);
    tau_x_hat_all = tau_x_mag_all;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% then move on to the exp reconstruction %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for (ite_pe = 1:max_pe_ite)
	

		    % give fixed magnitude measurements
		    % should i let it run until convergence

            %%% iteration in the measurement block
            tau_p_meas_2 = A_2_echo_nw.multSq(tau_x_hat_all);
            p_hat_meas_2 = A_2_echo_nw.mult(x_hat_all) - tau_p_meas_2 * s_hat_meas_2;


            for (ite_pe_est = 1:max_pe_est_ite)
                tau_w_2 = output_parameter_est(y, tau_w_2, p_hat_meas_2, tau_p_meas_2, kappa);
            end

            tau_s_meas_2 = 1 / (tau_w_2 + tau_p_meas_2);
            s_hat_meas_2 = (y - p_hat_meas_2) * tau_s_meas_2;

            tau_r_meas_2 = 1 / A_2_echo_nw.multSqTr(tau_s_meas_2);
            r_hat_meas_2 = x_hat_all + tau_r_meas_2 * A_2_echo_nw.multTr(s_hat_meas_2);

            tau_x_mag_meas_2 = tau_r_meas_2;
            x_mag_meas_2_complex = r_hat_meas_2;

            x_mag_meas_2 = abs(r_hat_meas_2);
            x_mag_meas_2(x_mag_meas_2<x_mag_min) = x_mag_min;

            tau_x_mag_meas = (tau_x_mag_meas_2*tau_x_mag_meas_1) / (tau_x_mag_meas_2 + tau_x_mag_meas_1);
            x_mag_meas_complex = (tau_x_mag_meas_2*x_mag_meas_1_complex + tau_x_mag_meas_1*x_mag_meas_2_complex) / (tau_x_mag_meas_2 + tau_x_mag_meas_1);
            x_mag_meas = abs(x_mag_meas_complex);
            x_mag_meas(x_mag_meas<x_mag_min) = x_mag_min;


            %%% iteration in the exp block
            tau_s_exp = 1 / (tau_x_mag_meas + tau_p_exp);
            s_hat_exp = (x_mag_meas - p_hat_exp) * tau_s_exp;

            tau_r_exp = 1 / A_x0_r2star.multSqTr(tau_s_exp);
            r_hat_exp = x0_psi + tau_r_exp * A_x0_r2star.multTr(s_hat_exp);

            % the input parameter estimation must be abs of r_hat_exp
            for (ite_pe_est = 1:max_pe_est_ite)
                lambda_x0_psi = input_parameter_est(abs(r_hat_exp), tau_r_exp, lambda_x0_psi, kappa);
            end

            x0_psi_pre = x0_psi;
            [x0_psi, tau_x0_psi] = input_function(r_hat_exp, tau_r_exp, lambda_x0_psi);

            % x0 should be positive always, enforce the lower bounds
            x0_tmp = A_wav_single.mult(x0_psi);
            x0_tmp = min(x0_tmp, x0_max);
            x0_tmp = max(x0_tmp, x0_min);
            x0_psi = A_wav_single.multTr(x0_tmp);

            cvg_x0_psi = norm(x0_psi(:)-x0_psi_pre(:), 'fro')/norm(x0_psi(:), 'fro');

            x0 = A_wav_single.mult(x0_psi);
            x0 = min(x0, x0_max);
            x0 = max(x0, x0_min);

            r2star = exp_lsq_fit_r2_componentwise_dict_2d(x_mag_meas, x0, exp_dict, exp_dict_sq_norm, r2_dict_val);
            r2star = min(r2star, r2star_max);
            r2star = max(r2star, r2star_min);
          
            set(A_x0_r2star, 'R2star', r2star);
            A_x0_r2star.estFrob();

            tau_p_exp = A_x0_r2star.multSq(tau_x0_psi);
            p_hat_exp = A_x0_r2star.mult(x0_psi) - tau_p_exp * s_hat_exp;

            %%% compute tau_x_hat_all and x_hat_all
            x_mag_exp = p_hat_exp;
            x_mag_exp = min(x_mag_exp, x_mag_max);
            x_mag_exp = max(x_mag_exp, x_mag_min);
            tau_x_mag_exp = tau_p_exp;

            tau_x_mag_all = (tau_x_mag_exp*tau_x_mag_meas) / (tau_x_mag_exp+tau_x_mag_meas);
            x_mag_all = (tau_x_mag_exp*x_mag_meas+tau_x_mag_meas*x_mag_exp) / (tau_x_mag_exp+tau_x_mag_meas);

            x_hat_all_pre = x_hat_all;
            x_hat_all = x_mag_all .* sign(x_mag_meas_complex);
            x_hat_all = x_hat_all_pre + damp_rate*(x_hat_all-x_hat_all_pre);
            tau_x_hat_all = tau_x_mag_all;

            cvg_x_hat_all = norm(x_hat_all(:)-x_hat_all_pre(:), 'fro')/norm(x_hat_all(:), 'fro');

		    cvg_gamp = max([cvg_x0_psi cvg_x_hat_all]);
		    if (cvg_gamp<cvg_thd) 
		    	break;
		    end
		    


        fprintf('Ite %d CVG PE: %d\t%d\n', ite_pe, cvg_x0_psi, cvg_x_hat_all)

        cvg_pe = max([ cvg_x0_psi cvg_x_hat_all ]);;
        if ((cvg_pe<cvg_thd)&&(ite_pe>2))
            break;
        end
	end
	
	res.x0_psi = x0_psi;
	res.tau_x0_psi = tau_x0_psi;
    res.r2star = r2star;

    res.x_hat_all = x_hat_all;
    res.tau_x_hat_all = tau_x_hat_all;
    res.s_hat_meas_2 = s_hat_meas_2;

    res.x_mag_exp = x_mag_exp;
    res.tau_x_mag_exp = tau_x_mag_exp;

    res.tau_p_exp = tau_p_exp;
    res.p_hat_exp = p_hat_exp;

    input_par.lambda_x0_psi = lambda_x0_psi; % Bernoulli parameter
    input_par.lambda_x_hat_psi = lambda_x_hat_psi;
    output_par.tau_w_1 = tau_w_1;
    output_par.tau_w_2 = tau_w_2;

end

function [x0_hat, tau_x0] = input_function(r_hat, tau_r, lambda)

    thresh = lambda * tau_r;
    x0_hat = max(0, abs(r_hat)-thresh) .* sign(r_hat);

    tau_x0 = tau_r;
    %tau_x0(abs(x0_hat)==0) = 0;
    tau_x0 = tau_r * length(x0_hat(abs(x0_hat)>0)) / length(x0_hat);

end


function lambda = input_parameter_est(r_hat, tau_r, lambda, kappa)

    lambda_pre = lambda;

    % do we need some kind of normalization here?
    dim_smp=length(r_hat);
    num_cluster=1;

    block_mat = zeros(dim_smp, num_cluster);
    for (i=1:num_cluster)
        block_mat(:,i) = 0.5/tau_r * (tau_r*lambda(i)-r_hat).^2;
    end

    block_mat_min = [];
    if (num_cluster==1)
        block_mat_min = block_mat;
    else
        block_mat_min = max(block_mat')';
    end

    block_mat = block_mat - block_mat_min; % subtract the minimum value of each row
    block_mat_two_exp = exp(block_mat);
    block_mat_one_exp = exp(-block_mat_min);

    block_mat_erfc = zeros(dim_smp, num_cluster);
    for (i=1:num_cluster)
        block_mat_erfc(:,i) = erfc(sqrt(0.5/tau_r)*(tau_r*lambda(i)-r_hat));
    end

    lambda_tmp_mat_0 = zeros(dim_smp, num_cluster);
    for (i = 1:num_cluster)
        lambda_tmp_mat_0(:,i) = lambda(i)/2 * block_mat_two_exp(:,i) .* block_mat_erfc(:,i);
    end

    % compute lambda
    % compute the first order derivative

    der_block = zeros(dim_smp, num_cluster);
    for (i=1:num_cluster)
        der_block(:,i) = ( sqrt(2*tau_r/pi) ./ erfcx(sqrt(0.5/tau_r)*(tau_r*lambda(i)-r_hat)) + r_hat - tau_r* lambda(i) );
    end

    fst_der_lambda = zeros(dim_smp, num_cluster);
    for (i = 1:num_cluster)
        fst_der_lambda(:,i) = 1/lambda(i) - der_block(:,i);
    end

    scd_der_lambda = zeros(dim_smp, num_cluster);
    for (i = 1:num_cluster)
        scd_der_lambda(:,i) = -1/lambda(i)/lambda(i) + ( tau_r + (r_hat - tau_r*lambda(i)).*der_block(:,i) ) - der_block(:,i).^2;
    end

    lambda_tmp_mat = lambda_tmp_mat_0;
    lambda_tmp_mat_sum = sum(lambda_tmp_mat, 2);
    for (i=1:num_cluster)
        lambda_tmp_mat(:,i) = lambda_tmp_mat(:,i) ./ ( lambda_tmp_mat_sum + eps);
    end

    lambda_tmp_mat_1 = lambda_tmp_mat;
    lambda_tmp_mat_2 = lambda_tmp_mat;
    %for (i=1:num_cluster)
    %    lambda_tmp_mat_1(:,i) = lambda_tmp_mat_1(:,i) .* (fst_der_lambda(:,i) - scd_der_lambda(:,i)*          lambda(i));
    %    lambda_tmp_mat_2(:,i) = lambda_tmp_mat_2(:,i) .* scd_der_lambda(:,i);
    %end
    %lambda_new = - sum(lambda_tmp_mat_1) ./ (sum(lambda_tmp_mat_2) + eps);   % to avoid division by 0
    %lambda_new = lambda_new';

    lambda_new = [];
    for (i=1:num_cluster)
        lambda_tmp_mat_1_tmp = lambda_tmp_mat_1(:,i) .* fst_der_lambda(:,i);
        lambda_tmp_mat_2_tmp = lambda_tmp_mat_2(:,i) .* scd_der_lambda(:,i);
        lambda_new_tmp = 0;
        if (sum(lambda_tmp_mat_2_tmp)<0)
            lambda_new_tmp = lambda(i) - sum(lambda_tmp_mat_1_tmp)/sum(lambda_tmp_mat_2_tmp);
        else
            if (sum(lambda_tmp_mat_1_tmp)>0)
                lambda_new_tmp = lambda(i)*1.1;
            else
                lambda_new_tmp = lambda(i)*0.9;
            end
        end
        lambda_new = [lambda_new; lambda_new_tmp];
    end

    lambda_new = max(lambda_new, 1e-12);  % necessary to avoid 0 which leads to NaN

    % lambda_new could be negative, causing problem
    lambda = lambda + kappa * (lambda_new - lambda);
end


function tau_w = output_parameter_est(y, tau_w, p_hat, tau_p, kappa)

    %tau_w_new = max(mean( abs(p_hat(:)-y(:)).^2 - tau_p(:) ), eps);   % make sure the variance is non-negative
    tau_w_new = mean(abs(p_hat(:)-y(:)).^2) + tau_p;
    %tau_w_new = mean( (p_hat-y).^2 - tau_p );
    %if (tau_w_new<0)
    %    tau_w_new = tau_w;
    %end

    tau_w = tau_w + kappa * (tau_w_new - tau_w);
end

