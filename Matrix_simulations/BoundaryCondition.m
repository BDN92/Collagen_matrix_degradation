% function to set Boundary conditions for spheres when they do diffusive motion

function [new_pos_s] = BoundaryCondition(new_pos_s, H, W, L)

% periodic BC in X direction
if new_pos_s(1) < 0
   new_pos_s(1) = L + new_pos_s(1); 
end
if new_pos_s(1) > L
   new_pos_s(1) = new_pos_s(1) - L;
end

% periodic BC in Y direction
if new_pos_s(2) < 0
   new_pos_s(2) = W + new_pos_s(2); 
end
if new_pos_s(2) > W
   new_pos_s(2) = new_pos_s(2) - W;
end

% periodic BC in Z direction
if new_pos_s(3) < 0
   new_pos_s(3) = H + new_pos_s(3); 
end
if new_pos_s(3) > H
   new_pos_s(3) = new_pos_s(3) - H;
end

end