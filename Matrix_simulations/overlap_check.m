% This overlap check function is related to Browdian dynamics simulation
% After each time step, when sphere copordinates are updated, then their
% overlap with fibers and other spheres are checked
% Cichocki-Hinsen algorith to simulate hindered diffusion is based on this overlap function 

function [condition] = overlap_check(id, new_pos_s, Nodes_Fibers, Fibers, dia_fiber, dia_sphere, sc, no_fibers)

condition = [];

xs = new_pos_s(1);
ys = new_pos_s(2);
zs = new_pos_s(3);

% check with fibers
for j = 1:no_fibers          
            
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

% check wih other spheres
for j  = 1:length(sc)
    if id ~= j
        dis = ((xs - sc(j,1))^2 + (ys - sc(j,2))^2 + (zs - sc(j,3))^2)^0.5;
            if dis > dia_sphere
                condition = [condition; 0];
            else
                condition = [condition; 1];
            end
     end
end

end