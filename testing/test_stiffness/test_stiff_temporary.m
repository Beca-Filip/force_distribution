ttt = 10;
sss = 2;
lll = 2;

figure;
plot(squeeze(data.K(:,:,:,ttt, sss, lll)).')
plot(squeeze(data.Kref(:,:,:,ttt, sss, lll)).', '--', 'linewidth', 2)