function r2star = exp_lsq_fit_3d(x0, x_mag_meas, echo_time)

    x0(x0<eps) = eps;
    x_mag_meas(x_mag_meas<eps) = eps;

    x0_log = log(x0);
    x_mag_meas_log = log(x_mag_meas);
    x_mag_meas_sq = x_mag_meas.^2;

    echo_time_sq = echo_time.^2;
    
    r2star_num = zeros(size(x0));
    r2star_den = zeros(size(x0));
    for (i=1:length(echo_time))
        r2star_num = r2star_num + echo_time(i)*x_mag_meas_sq(:,:,:,i).*(x0_log-x_mag_meas_log(:,:,:,i));
        r2star_den = r2star_den + echo_time_sq(i)*x_mag_meas_sq(:,:,:,i);
    end

    r2star_den(r2star_den<eps) = eps;
    r2star = r2star_num ./ r2star_den;
end
