% perform reconstruction in the fourier domain
% start with the wavelet implementation to make sure this framework will work, testing 

load('./complex_image.mat')
load('./sensitivity_map_3d.mat')

sx = size(complex_image,1);
sy = size(complex_image,2);
sz = size(complex_image,3);
Ne = size(complex_image,4);
Nc = size(maps_3d,4);

for (i=1:Ne)
    data = zeros(sx,sy,sz,Nc);
    for (j=1:Nc)
        data_tmp = fftshift(fftn( complex_image(:,:,:,i) .* maps_3d(:,:,:,j) ));
        data(:,:,:,j) = data_tmp;
    end

    data = single(data);

    save(strcat('./echo_', num2str(i), '_3d.mat'), 'data', '-v7.3')
end
