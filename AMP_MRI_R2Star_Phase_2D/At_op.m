function X = At_op(y, s_vect, mat_sz) 
    X_tmp = zeros(mat_sz);
    X_tmp(s_vect) = y;
    X = ifft2(ifftshift(X_tmp)); 
end

