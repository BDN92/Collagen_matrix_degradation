clear
clc

% Simulating Brownian motion of MMPs through 
% collagen matix of uniform fibers, i.e., they are of same length and diameter

% simulation box 
% all length dimensions are scaled by 1 micron
H = 5;    
W = 5;     
L = 5;    

x = [0 L L 0];                                                        
y = [0 W];                                                                  
z = [0 0 H H];

fibril_fraction = 6 * 10^(-3);                % setting fibril fraction
vol_fibrils = (H * W * L) * fibril_fraction;

length_fiber = 2;
dia_fiber = 0.1;
no_fibers = vol_fibrils / (pi/4 * dia_fiber^2 * length_fiber);
no_fibers = round(no_fibers);

No_spheres = 1500;  % number of enzymes
dia_sphere = 0.01;  % hydrodynamic diameter of enzyme in micron

% generate fibrils and spheres
[Nodes_Fibers, Fibers, spheres_coordinates] = many_fibrils_gen(H, W, L, no_fibers, length_fiber, dia_fiber, No_spheres, dia_sphere);

%%%%%% BD simulation %%%%%%%
No_timesteps = 10^6;
tau = 10^(-6);   % (in second) simulation time step
D = 74;        % (scaled by 1 micron^2/s) diffusivity bacterial collagenase 74 * 10^(-12) m^2/s (Schulz & Anseth, Soft Matter, 2013)

sc_initial = spheres_coordinates;
sc = spheres_coordinates;
sc_displacement = zeros(No_spheres,3); 

del_t_sampling = 50;
spheres_coordinates_time = cell(No_timesteps/del_t_sampling,2);
spheres_displacement_time = cell(No_timesteps/del_t_sampling,2);

for i = 1:No_timesteps

    k = sqrt(2 * 1 * D * tau);     % increment factor in one dimension

    % updating position
    for j = 1:No_spheres
        id = j;
        old_pos_s = sc(id,:);
        ds = k * randn(1,3);
        new_pos_s = old_pos_s + ds;

        % periodic BC in y and z directions and reflective BC in x direction
        [new_pos_s] = BoundaryCondition(new_pos_s, H, W, L);

        % checking new_pos_s overlapping with fibers and other spheres
        [condition] = overlap_check(id, new_pos_s, Nodes_Fibers, Fibers, dia_fiber, dia_sphere, sc, no_fibers);

        if mean(condition) == 0   % no overlap and update position
            sc(id,:) = new_pos_s;
            sc_displacement(id,:) = ds;
        end
    end

    % data sampling
    if mod(i,del_t_sampling) == 0
        spheres_coordinates_time{i/del_t_sampling,1} = i;
        spheres_coordinates_time{i/del_t_sampling,2} = sc;
        spheres_displacement_time{i/del_t_sampling,1} = i;
        spheres_displacement_time{i/del_t_sampling,2} = sc_displacement;

        % tracking how many time steps are elapsed
        disp(['No. time steps elapsed:' num2str(i)]);
    end

    if mod(i,50000) == 0
        outputname = ['lf' num2str(length_fiber) '_df' num2str(dia_fiber) '_F' num2str(no_fibers) '_E' num2str(No_spheres) '_dE' num2str(dia_sphere) '_' num2str(i) '.mat'];        
        save(outputname, 'H', 'W', 'L', 'no_fibers', 'length_fiber', 'dia_fiber', 'No_spheres', 'dia_sphere', 'Nodes_Fibers', 'Fibers', 'tau', 'D', 'spheres_coordinates_time', 'spheres_displacement_time' )
    
    end
end

% plotting spheres and fibers
 spheres_coordinates_final = spheres_coordinates_time{No_timesteps/del_t_sampling,2};
 Plot_Fiber(x,y,z,Nodes_Fibers,Fibers,dia_fiber);
 hold on
 Plot_spheres(No_spheres, dia_sphere, spheres_coordinates)

