% function file to plot solid spheres

function Plot_spheres(No_spheres, dia_sphere, spheres_coordinates)

for i = 1:No_spheres
    [p q r] = sphere;
    s = surf(p * dia_sphere + spheres_coordinates(i,1), q * dia_sphere + spheres_coordinates(i,2), r * dia_sphere + spheres_coordinates(i,3));
    set(s,'FaceColor','r','EdgeColor','none','AmbientStrength',.5);
    view(60,30)
end

end