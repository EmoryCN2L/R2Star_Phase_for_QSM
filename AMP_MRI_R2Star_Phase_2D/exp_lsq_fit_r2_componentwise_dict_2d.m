 function V=exp_lsq_fit_r2_componentwise_dict_2d(Et, H_exp, exp_dict, exp_dict_sq_norm, r2_dict_val)
 
     % solve the least square problem for R2 only while keeping H_exp fixed
 
     % Et        : the recovered magnitude echo images of size sx by sy by Ne
     %           : sx is the resolution in the x direction
     %           : sy is the resolution in the y direction
     %           : Ne is the number of echoes
 
     % par       : various parameters
 
     % Assigning parameres according to par
 
     % maximum and minimum value of Rt_abs
 
     %if (maxiter<2000)
     %    fprintf('Increasing maxiter to 2000!\n')
     %    maxiter = 2000;
     %end
 
     sx = size(Et,1);
     sy = size(Et,2);
     echo_num = size(Et,3);
 
     H_exp_min = eps;
     Et_min = eps;
     V_min = 0;  % or eps???
 
     % initialize H_exp and V
     Et = max(Et, Et_min);
     
     % initialize H, H_exp, and V
     H_exp = max(H_exp, H_exp_min);
 
    Et_process = Et;
    for (i=1:size(Et,3))     
        Et_process(:,:,i) = Et_process(:,:,i) ./ H_exp;
    end
    
    Et_process_sq_norm = sum(Et_process.^2,3);
    
    min_distance_to_dict = zeros(sx,sy);
    V = zeros(sx,sy);
    for (i=1:size(exp_dict,2))
        %if (mod(i,500)==0)
        %fprintf('%d\n', i)
        %end
        exp_dict_seq = exp_dict(:,i);
        distance_to_dict = zeros(sx,sy);
        for (j=1:length(exp_dict_seq))
            distance_to_dict = distance_to_dict + Et_process(:,:,j)*exp_dict_seq(j);
        end
        distance_to_dict = Et_process_sq_norm+exp_dict_sq_norm(i)-2*distance_to_dict;
        
        if (i==1)
            V = repmat(r2_dict_val(i),[sx sy]);
            min_distance_to_dict = distance_to_dict;
        else
            V(min_distance_to_dict>distance_to_dict) = r2_dict_val(i);
            min_distance_to_dict = min(min_distance_to_dict, distance_to_dict);
        end            
    end

end
