eigThresh_im = 0;
Nc=32;
for (i=1:208)
    fprintf('%d\n',i)
    load(strcat('./result_sen_map_3d/sen_map_', num2str(i), '_W.mat'))
    load(strcat('./result_sen_map_3d/sen_map_', num2str(i), '_M.mat'))
    M=squeeze(M);
    W=squeeze(W);
    maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_im,[1,1,Nc]);
    for (j=1:Nc)
        maps(:,:,j) = fftshift(squeeze(maps(:,:,j)));
    end 
    save(strcat('./result_sen_map_3d/sen_map_', num2str(i), '.mat'), 'maps')
end

