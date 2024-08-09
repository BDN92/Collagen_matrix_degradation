% function generating non-overlapping fibers and spheres
% using fiber generatgion and sphere generation functions

function [Nodes_Fibers, Fibers, spheres_coordinates] = many_fibrils_gen(H, W, L, no_fibers, length_fiber, dia_fiber, No_spheres, dia_sphere);


Orientation = []; % empty orientation represents random direction
Ndiv = 1;
x = [0 L L 0];                                                        
y = [0 W];                                                                  
z = [0 0 H H];


[Nodes_Fibers, Fibers] = Generate_Fiber(x,y,z,length_fiber,no_fibers,dia_fiber,Orientation,Ndiv); 

[spheres_coordinates] = Generate_Sphere(x,y,z, No_spheres, dia_sphere, dia_fiber, Nodes_Fibers, Fibers);

end