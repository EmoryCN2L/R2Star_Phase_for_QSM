function X = At_op_3d_cylinder(y, s_vect, mat_sz) 
    X_tmp = zeros(mat_sz);
    for (i=1:mat_sz(1))
        X_tmp_tmp = zeros(mat_sz(2), mat_sz(3));
        X_tmp_tmp(s_vect) = y(i,:);
        X_tmp(i,:,:) = X_tmp_tmp;
    end
    X = ifftn(ifftshift(X_tmp));    
end

