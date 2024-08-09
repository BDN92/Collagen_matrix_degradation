%%%% checking whether generated fibers overlap or not
%%%% using a logic based approach

div = 100;
x_points = zeros(no_fibers,div);
y_points = zeros(no_fibers,div);
z_points = zeros(no_fibers,div);

for i = 1:no_fibers
   
        x1 = Nodes_Fibers(Fibers(i,1),1);
        y1 = Nodes_Fibers(Fibers(i,1),2);
        z1 = Nodes_Fibers(Fibers(i,1),3);

        x2 = Nodes_Fibers(Fibers(i,2),1);
        y2 = Nodes_Fibers(Fibers(i,2),2);
        z2 = Nodes_Fibers(Fibers(i,2),3);

        xq = linspace(x1,x2,div);
        yq = linspace(y1,y2,div);
        zq = linspace(z1,z2,div);

        x_points(i,:) = xq;
        y_points(i,:) = yq;
        z_points(i,:) = zq;
end

logical = ones(no_fibers,no_fibers);

for i = 1:no_fibers 
    for j = 1:no_fibers
        if i ~= j
            xx1 = x_points(i,:)';
            xx2 = x_points(j,:)';

            yy1 = y_points(i,:)';
            yy2 = y_points(j,:)';

            zz1 = z_points(i,:)';
            zz2 = z_points(j,:)';

            dis = zeros(div,div);

            for k = 1:div
                for m = 1:div
                    d = ((xx1(k,1) - xx2(m,1))^2 + (yy1(k,1) - yy2(m,1))^2 + (zz1(k,1) - zz2(m,1))^2)^0.5; 
                    dis(k,m) = d;
                end
            end

            if min(min(d)) > dia_fiber
               logical(i,j) =  0;
            end    
        end  
    end         
end

overlap = 0;
for i = 1:length(logical)
    for j = 1:length(logical)
        if logical(i,j) == 1 && i~=j
            overlap = overlap+1;
        end
    end
end

if overlap == 0
    overlap
    disp("No Overlap Among all fibers")
else
    overlap
    disp("There are overlaps")
end






