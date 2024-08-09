
% condition function file to check overlap
function doOverlap = box_overlap(x1,y1,z1,x2,y2,z2)

  doOverlap = false;
  
  flagx = false; 
  flagy = false; 
  flagz = false;

  x_min = x1(1); x_max = x1(2); 
  x_min2 = x2(1); x_max2 = x2(2);
  
  y_min = y1(1); y_max = y1(2); 
  y_min2 = y2(1); y_max2 = y2(2);
  
  z_min=z1(1); z_max=z1(2); 
  z_min2=z2(1); z_max2=z2(2);

  if (x_min <= x_max2) && (x_min2 <= x_max)
     flagx=true;
  end
  
  if (y_min <= y_max2) && (y_min2 <= y_max)
     flagy = true;
  end
  
  if (z_min <= z_max2) && (z_min2 <= z_max)
     flagz = true;
  end

  if flagx == true && flagy == true && flagz == true
     doOverlap = true;
  end

end