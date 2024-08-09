
% function file to generate mesh points required to check overlap
function C = Cuboid(x,y,z,N_div_x,N_div_y,N_div_z,xb,yb,zb)

x1 = min(x); x2=max(x);
y1 = min(y); y2=max(y);
z1 = min(z); z2=max(z);

dx = (x2-x1)/N_div_x;
dy = (y2-y1)/N_div_y;
dz = (z2-z1)/N_div_z;

xgv = floor((xb(1)-x1)/dx):1:floor((xb(2)-x1)/dx);
ygv = floor((yb(1)-y1)/dy):1:floor((yb(2)-y1)/dy);
zgv = floor((zb(1)-z1)/dz):1:floor((zb(2)-z1)/dz);

[Y,X,Z] = meshgrid(ygv,xgv,zgv);

V = X + Y * N_div_x + Z * N_div_x * N_div_y + 1;

C = reshape(V,1,numel(V));

end