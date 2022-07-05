ttt = 10;
sss = 3;
lll = 1;

figure;
plot(squeeze(data.K(:,:,:,ttt, sss, lll)).')
plot(squeeze(data_copy.K(:,:,:,ttt, sss, lll)).', '--', 'linewidth', 2)