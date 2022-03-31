% construct linear measurement operator and its transpose operator
function y = A_op(X, s_vect)
    X_fft = fftshift(fft2(X)); 
    y=X_fft(s_vect);
end

