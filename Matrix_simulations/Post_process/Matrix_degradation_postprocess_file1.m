clear
clc

% post processing to obtain enzyme per fiber and to generate enzyme distribution
% using enzyme values, predict diameters
% enzyme per fiber is obtained by minimum distance of enzyme from a fiber

% first check steady state
% for that run the script 'check_find_steady_state.m'
% if steady state interval is found
% then run the next lines after choosing appropriate time interval


% below is an example file
% for a matrix having uniform fibrils
% generate this example file below by running Main_collagen.m
load('lf2_df0.1_F48_E1500_dE0.01_300000.mat')    % load file generated from Brownian Dynamics

%% time sampling and time averaging 
sampling = 200;   % number of sample to be averaged 
sg = 6000;        % total number of samples generated from time t = 0 
post_process_sc = cell(sampling,2);
for i = 1:sampling
    post_process_sc(i,1) = spheres_coordinates_time(sg - sampling + i - 1, 1);     
    post_process_sc(i,2) = spheres_coordinates_time(sg - sampling + i - 1, 2);
end

fiber_enzyme_stat = zeros(no_fibers,sampling);

for i = 1:sampling
    sc = post_process_sc{i,2};

    for j = 1:No_spheres
        dis = zeros(no_fibers,1);
        for k = 1:no_fibers

            l = Nodes_Fibers(Fibers(k,2),1) - Nodes_Fibers(Fibers(k,1),1);
            m = Nodes_Fibers(Fibers(k,2),2) - Nodes_Fibers(Fibers(k,1),2);
            n = Nodes_Fibers(Fibers(k,2),3) - Nodes_Fibers(Fibers(k,1),3);
        
            mod_unitvec = (l^2 + m^2 +n^2)^0.5;
        
            xfc = (Nodes_Fibers(Fibers(k,2),1) + Nodes_Fibers(Fibers(k,1),1))/2;
            yfc = (Nodes_Fibers(Fibers(k,2),2) + Nodes_Fibers(Fibers(k,1),2))/2;
            zfc = (Nodes_Fibers(Fibers(k,2),3) + Nodes_Fibers(Fibers(k,1),3))/2;

            xs = sc(j,1);
            ys = sc(j,3);
            zs = sc(j,2);
            xdis = xs - Nodes_Fibers(Fibers(k,1),1);
            ydis = ys - Nodes_Fibers(Fibers(k,1),2);
            zdis = zs - Nodes_Fibers(Fibers(k,1),3);

            det = ((ydis * n - zdis * m)^2 + (xdis * n - zdis * l)^2 + (xdis * m - ydis * l)^2)^0.5;
            
            % min perpendicular distance 
            distance1 = det/mod_unitvec;

            % center-to-center distance
            distance2 = ((xs - xfc)^2 + (ys - yfc)^2 + (zs - zfc)^2)^0.5;

            % conditions to choose distance --- case 1
            base = (distance2^2 - distance1^2)^0.5;
            if base < length_fiber/2
               dis(k,1) = distance1 - dia_fiber/2;
            end
            % case 2
            if base > length_fiber/2
               base2 = (distance2^2 - distance1^2)^0.5 - length_fiber/2;
               case2_dis = (base2^2 + (distance1 - dia_fiber/2)^2)^0.5;
               dis(k,1) = case2_dis;
            end

        end
        min_dis = min(dis);
        [row, col] = find(dis==min_dis);
        fiber_enzyme_stat(row,i) = fiber_enzyme_stat(row,i) + 1; 
    end
end
fiber_enzyme_av_dist = zeros(no_fibers,1);      % enzyme per fibril
fiber_enzyme_avstd_dist = zeros(no_fibers,1);    % error 
for i = 1:no_fibers
    fiber_enzyme_av_dist(i,1) = mean(fiber_enzyme_stat(i,:));
    fiber_enzyme_avstd_dist(i,1) = std(fiber_enzyme_stat(i,:));
end
fiber_enzyme_avstd_dist = 1.984 * fiber_enzyme_avstd_dist/200^0.5;
fiberid = (1:1:no_fibers)';

% % histogram fit function
% h1=histogram(fiber_enzyme_av_dist,20, 'Normalization', 'pdf','FaceColor', 'b');
% h1.BinCounts(h1.BinCounts>0) = h1.BinCounts(h1.BinCounts>0);
% x = h1.BinEdges(1:end-1)+h1.BinWidth/2;
% y = h1.Values;
% % Lognormal function
% %lognormal = @(mu, sig, scale, x) scale./(x*sig*sqrt(2*pi)) .* exp(-((log(x)-mu).^2)./(2*sig^2));
% lognormal = @(mu, sig, x) 1./(x*sig*sqrt(2*pi)) .* exp(-((log(x)-mu).^2)./(2*sig^2)); 
% %x0 = [mean(x), range(x), sum(h1.Values*h1.BinWidth)];
% x0 = [mean(log(x)), std(log(x))];
% % Fitting with a different optimization algorithm
% f1 = fit((x(:)), y(:), lognormal, 'StartPoint', x0);
% hold on;
% p = plot(f1);
% p.Color = 'b';
% p.LineWidth = 2;
% hold on;


%% Run single fibril model using enzyme values, 
% and predict diameters of fibrils
df_0 = dia_fiber * 1000;     % converting from micrometer to nanometer 
lf_0 = length_fiber * 1000;

DF_final = zeros(no_fibers,1);   

for i = 1:no_fibers
    NE_0 = fiber_enzyme_av_dist(i,1);
    time = 3600;
    [t, df_final] = single_fibril_thickness_reduction_kinetics(df_0, lf_0, NE_0, time);
    DF_final(i,1) = df_final(length(t),1);
end
DF = DF_final/1000;  % converting unit back to micrometer

%%%% total mass lost estimation
initial_mass = no_fibers * (pi/4) * dia_fiber^2 * length_fiber;
final_mass = 0;
for i = 1:no_fibers
    final_mass = final_mass + (pi/4) * (DF_final(i,1)/1000)^2 * length_fiber;
end

result_data = [fiberid fiber_enzyme_av_dist fiber_enzyme_avstd_dist DF];
mass_lost = (initial_mass - final_mass)/initial_mass;


%% Note: If one fibril is completely dissolved before set time, 
% then enzymes corresponding to that fibrils must be assigned to
% other fibrils based on next nearest distance
% this is not implemented in this current code
% the current code is valid for the cases where fibrils do not dissolve
% completely before set time, thus reassignment of enzymes does not require
% if fibrils dissolve, reassignment of enzymes is required








