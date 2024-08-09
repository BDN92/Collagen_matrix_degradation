% function file generating nodes of fibers and fiber ids

function [Nodes_Fibers, Fibers]=Generate_Fiber(x,y,z,length_fiber,no_fibers,dia_fiber,Orientation,Ndiv)

% disp('Generate Fiber');
% disp('-----------------------');

% input

% x,y,z: as vectors defining the simulation box dimension 
% L: Length of fiber
% Dfiber = Diameter of fiber
% N: Number of fibers
% Orientation: Orientation of fibers as [l; m; n] column vector where l, m, and n are direction cosines of orientation in x, y, and z directions, respectively
% Ndiv: Number of fiber mesh divisions

N = no_fibers;
L = length_fiber;
DFiber = dia_fiber;

N_el = N;

lx = max(x)-min(x);   % length of sample in x-axis
ly = max(y)-min(y);   % length of sample in y-axis
lz = max(z)-min(z);   % length of sample in z-axis

N_div_x = round((N_el/(lx*ly*lz))^(1/3)*lx);   % Number of portions in x-axis
N_div_y = round((N_el/(lx*ly*lz))^(1/3)*ly);   % Number of portions in y-axis
N_div_z = round((N_el/(lx*ly*lz))^(1/3)*lz);   % Number of portions in z-axis

C_max = N_div_x * N_div_y * N_div_z;             % Number of small cuboids (inside simulation box) to insert fibers 

for l = 1:1:C_max
    CuboidsF{l}=[];
end

xxF = zeros(N,2);
yyF = zeros(N,2);
zzF = zeros(N,2);

XBF = zeros(N,2);
YBF = zeros(N,2);
ZBF = zeros(N,2);

Nodes_Fibers = [];   
Fibers = [];

tic;

% start inserting fibers one by one
for i = 1:1:N

   if toc>2        
      disp([num2str(i/N*100) '% Completed'])
      tic
   end


while 1

x1 = min(x)+(max(x)-min(x))*rand(1);
y1 = min(y)+(max(y)-min(y))*rand(1);
z1 = min(z)+(max(z)-min(z))*rand(1);

if isempty(Orientation)
    theta = 2*pi*rand(1);
    v = 2*rand(1)-1;
    c2 = [0 0 0]';
    c2(1) = x1 + L*sqrt(1-v^2)*cos(theta);
    c2(2) = y1 + L*sqrt(1-v^2)*sin(theta);
    c2(3) = z1 + L*v;
else
    c2=[x1 y1 z1]' + L * Orientation/norm(Orientation)*sign(rand(1)*2-1);       
end

x2 = c2(1);
y2 = c2(2);
z2 = c2(3);

xq = linspace(x1,x2,10);
yq = linspace(y1,y2,10);
zq = linspace(z1,z2,10);

% Check if fibers overlap with boundaries    
if all(inpolygon(xq,zq,x,z)) && y2 <= y(2) && y2 >= y(1)  

% Find vectors perpendicular to fiber orientation 
Direction1(1,1) = x2-x1;
Direction1(2,1) = y2-y1;
Direction1(3,1) = z2-z1;
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

l3 = m1*n2 - m2*n1;
m3 = -l1*n2+l2*n1;
n3 = l1*m2-l2*m1;
Direction3 = [l3 m3 n3]; 
Direction3 = Direction3/norm(Direction3);

c = 0;
xx = zeros(218,1); yy = zeros(218,1); zz = zeros(218,1); 
for jk=1:1:12   
    Direction=Direction2*cos(2*pi*jk/12)-Direction3*sin(2*pi*jk/12);
    for k=1:1:10 
        c=c+1;    
        xx(c)=xq(k)+DFiber/2*Direction(1);
        yy(c)=yq(k)+DFiber/2*Direction(2);
        zz(c)=zq(k)+DFiber/2*Direction(3);  
    end
    for p=1:1:4
       c=c+1;
       xx(c)=xq(1)+DFiber/2*Direction(1)*p/5;
       yy(c)=yq(1)+DFiber/2*Direction(2)*p/5;
       zz(c)=zq(1)+DFiber/2*Direction(3)*p/5;  
       c=c+1;
       xx(c)=xq(10)+DFiber/2*Direction(1)*p/5;
       yy(c)=yq(10)+DFiber/2*Direction(2)*p/5;
       zz(c)=zq(10)+DFiber/2*Direction(3)*p/5;  
    end
end

c = c+1;
xx(c) = xq(1);
yy(c) = yq(1);
zz(c) = zq(1); 
c = c+1;
xx(c) = xq(10);
yy(c) = yq(10);
zz(c) = zq(10); 

%Check with fibers that lie within the same cuboid
xb = [min([x1 x2])-DFiber/2 max([x1 x2])+DFiber/2]; 
if xb(1)<min(x); xb(1)=min(x); end 
if xb(2)>max(x); xb(2)=max(x)-0.000001; end 

yb = [min([y1 y2])-DFiber/2 max([y1 y2])+DFiber/2]; 
if yb(1)<min(y); yb(1)=min(y); end
if yb(2)>max(y); yb(2)=max(y)-0.000001; end

zb = [min([z1 z2])-DFiber/2 max([z1 z2])+DFiber/2]; 
if zb(1)<min(z); zb(1)=min(z); end 
if zb(2)>max(z); zb(2)=max(z)-0.000001; end

doOverlap4 = false;
C = Cuboid(x,y,z,N_div_x,N_div_y,N_div_z,xb,yb,zb);
CC = []; %Vector of all fibers that lie in the same cuboids
for l=1:1:length(C)
    CC=union(CC,CuboidsF{C(l)});     
end
for l=1:1:length(CC)
    doOverlap2 = box_overlap(xb,yb,zb,XBF(CC(l),:),YBF(CC(l),:),ZBF(CC(l),:)); 
    if doOverlap2
       doOverlap3 = false;
       I = CheckPointCylinder([xxF(CC(l),1) ; yyF(CC(l),1) ; zzF(CC(l),1)],[xxF(CC(l),2) ; yyF(CC(l),2) ; zzF(CC(l),2)], DFiber, [xx yy zz]');
       if any(I == 1) 
          doOverlap3 = true;
       end
       if doOverlap3
          doOverlap4 = true;
          break
       end
    end
end
    
if not(doOverlap4)
   if not(all(inpolygon(xx,zz,x,z)) && all(yy<=y(2)) && all(yy>=y(1))) 
       doOverlap4 = true; 
   end  
end

if not(doOverlap4)
   break;
end

end    

end

    xxF(i,:)=[x1 x2];
    yyF(i,:)=[y1 y2];
    zzF(i,:)=[z1 z2];

    XBF(i,:)=xb;
    YBF(i,:)=yb;
    ZBF(i,:)=zb;

    for l=1:1:length(C)
         CuboidsF{C(l)}=[CuboidsF{C(l)} i];  
    end

    for j=1:1:Ndiv+1
        Nodes_Fibers((Ndiv+1)*(i-1)+j,1) = x1+(x2-x1)*(j-1)/Ndiv;
        Nodes_Fibers((Ndiv+1)*(i-1)+j,2) = y1+(y2-y1)*(j-1)/Ndiv;
        Nodes_Fibers((Ndiv+1)*(i-1)+j,3) = z1+(z2-z1)*(j-1)/Ndiv;    
    end

    for j=1:1:Ndiv
        Fibers(Ndiv*(i-1)+j,1)=(Ndiv+1)*(i-1)+j;   
        Fibers(Ndiv*(i-1)+j,2)=(Ndiv+1)*(i-1)+j+1;   
    end

end

disp('  ');

end
    