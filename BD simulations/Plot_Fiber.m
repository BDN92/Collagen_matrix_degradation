% function file to plot fibers 

function Plot_Fiber(x,y,z,Nodes_Fibers,Fibers,dia_fiber)

DFiber = dia_fiber;

figure;

% Plot simulation box
m = length(x);
xr = [x x(1)];
yr1 = y(1) * ones(m+1,1);
zr = [z z(1)];
plot3(xr,yr1,zr,'k','LineWidth',3);
hold on
yr2 = y(2) * ones(m+1,1);
plot3(xr,yr2,zr,'k','LineWidth',3);
hold on
for i = 1:1:m
    plot3([x(i) x(i)],y,[z(i) z(i)],'k','LineWidth',3);    
end



hold on


% Plot Fibers
for i = 1:1:size(Fibers,1)

x1 = Nodes_Fibers(Fibers(i,1),1);
y1 = Nodes_Fibers(Fibers(i,1),2);
z1 = Nodes_Fibers(Fibers(i,1),3);

x2 = Nodes_Fibers(Fibers(i,2),1);
y2 = Nodes_Fibers(Fibers(i,2),2);
z2 = Nodes_Fibers(Fibers(i,2),3);
  
Direction1(1,1) = x2 - x1;
Direction1(2,1) = y2 - y1;
Direction1(3,1) = z2 - z1;
Direction1 = Direction1/norm(Direction1);

l1 = Direction1(1,1); 
m1 = Direction1(2,1); 
n1 = Direction1(3,1);

if abs(n1) < 1e-6, n1 = 1e-6; end

l2 = rand(1); 
m2 = rand(1); 
n2 = (-l1*l2-m1*m2)/n1;
Direction2 = [l2 m2 n2]; 
Direction2 = Direction2/norm(Direction2);

l3 = m1 * n2 - m2 * n1;
m3 = -l1 * n2 + l2 * n1;
n3 = l1 * m2 - l2 * m1;
Direction3 = [l3 m3 n3]; 
Direction3 = Direction3/norm(Direction3);

xq = linspace(x1,x2,10);
yq = linspace(y1,y2,10);
zq = linspace(z1,z2,10);

c = 0;
xx = zeros(10,13); yy=zeros(10,13); zz=zeros(10,13); 

for jk = 1:1:13   
    Direction = Direction2*cos(2*pi*jk/12)-Direction3*sin(2*pi*jk/12);
    for k = 1:1:10 
        c = c+1;    
        xx(k,jk)=xq(k)+DFiber/2*Direction(1);
        yy(k,jk)=yq(k)+DFiber/2*Direction(2);
        zz(k,jk)=zq(k)+DFiber/2*Direction(3);  
     end
end

h=surface(xx, yy, zz);

set(h,'FaceColor',"#FF8800",'EdgeColor','none','AmbientStrength',.5);

end

hold('off');

daspect([1 1 1]);
camlight left;
lighting gouraud;
hold off
view(60,30)

end