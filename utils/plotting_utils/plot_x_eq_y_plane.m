function plot_x_eq_y_plane(xy_range, z_range)

xylo = min(xy_range);
xyhi = max(xy_range);

zlo = min(z_range);
zhi = max(z_range);

c1 = [xyhi*[1;1];zhi];
c2 = [xyhi*[1;1];zlo];
c3 = [xylo*[1;1];zlo];
c4 = [xylo*[1;1];zhi];
C = [c1, c2, c3, c4];

x = C(1, :).';
y = C(2, :).';
z = C(3, :).';

fill3(x,y,z,'b','FaceAlpha',0.3);

end