% function file generating coordinates of spheres

function [spheres_coordinates] = Generate_Sphere(x,y,z, No_spheres, dia_sphere, dia_fiber, Nodes_Fibers, Fibers)

% generating spheres at random positions non overlapping with fibers
spheres_coordinates = zeros(No_spheres,3);

factor_x = 1;
factor_y = 1;
factor_z = 1;

for i = 1:No_spheres

    while 1
        xs = (max(x)/factor_x - min(x)) * rand(1);
        ys = (max(y)/factor_y - min(y)) * rand(1);
        zs = (max(z)/factor_z - min(z)) * rand(1);
        
        condition = [];
        for j = 1:length(Fibers)
            l = Nodes_Fibers(Fibers(j,2),1) - Nodes_Fibers(Fibers(j,1),1);
            m = Nodes_Fibers(Fibers(j,2),2) - Nodes_Fibers(Fibers(j,1),2);
            n = Nodes_Fibers(Fibers(j,2),3) - Nodes_Fibers(Fibers(j,1),3);

            mod_unitvec = (l^2 + m^2 +n^2)^0.5;

            % finding distance between sphere and cylinder
            % line passing through (xs,ys,zs) and parallel to (l,m,n)
            % cross-product provides distance 
            xdis = xs - Nodes_Fibers(Fibers(j,1),1);
            ydis = ys - Nodes_Fibers(Fibers(j,1),2);
            zdis = zs - Nodes_Fibers(Fibers(j,1),3);

            det = ((ydis * n - zdis * m)^2 + (xdis * n - zdis * l)^2 + (xdis * m - ydis * l)^2)^0.5;
            distance = det/mod_unitvec; 

            if distance > (dia_sphere + dia_fiber)/2
                condition = [condition; 0];    % 0 for no overlap
            else 
                condition = [condition; 1];    % 1 for overlap
            end
        end

        for j  = 1:length(spheres_coordinates)
            dis = ((xs - spheres_coordinates(j,1))^2 + (ys - spheres_coordinates(j,2))^2 + (zs - spheres_coordinates(j,3))^2)^0.5;
            if dis > dia_sphere
                condition = [condition; 0];
            else
                condition = [condition; 1];
            end
        end


%         void =  isempty(spheres_coordinates);
%         if void == 0 
%             for j = 1:length(spheres_coordinates)
%                  dis = ((xs - spheres_coordinates(j,1))^2 + (ys - spheres_coordinates(j,2))^2 + (zs - spheres_coordinates(j,3))^2)^0.5;
%                   if dis > dia_sphere
%                      condition = [condition; 0];
%                   else
%                      condition = [condition; 1];
%                   end
%             end
%         end

        if mean(condition == 0)   % if all zero,  then no overlap with any fiber and spheres
            break;
        end

    end

    spheres_coordinates(i,:) = [xs ys zs];
end



end