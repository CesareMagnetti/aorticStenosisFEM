function [] = fem_surface_plot(var, field, varargin)

  n = 0;
  newfigure = 0;
  boundingbox=0;
  while(n < nargin-2)
    n = n + 1;
    if(strcmpi(varargin{n},'newfigure'))
      newfigure = 1;
    elseif(strcmpi(varargin{n},'BoundingBox'))
      n = n + 1; boundingbox=1; pts = varargin{n};
    end
  end

  if(newfigure)
    figure    
  end
  if(boundingbox)
    plot(pts(:,1),pts(:,2),'w.'); hold
  end
  
  if(strcmp(var.e.etype,'quadrilateral') == 1)
      if(var.e.p == 1)
          id = [1 2 3; 2 3 4];
      elseif(var.e.p == 2)
          id = [1 5 7; 1 6 7; 2 5 7; 2 7 8; 3 6 7; 3 9 7; 4 9 7; 4 8 7];         
      end
  elseif(strcmp(var.e.etype,'triangle') == 1)
      if(var.e.p == 1)
          id = [1 2 3];
      elseif(var.e.p == 2)
          id = [1 4 6; 1 5 6; 2 4 6; 3 5 6];         
      end      
  else
      error('Unknown element type')
  end
  
  for i = 1:size(id,1)
    trisurf(var.t(:,id(i,:)), var.x(:,1), var.x(:,2), field,'Facecolor','interp','LineStyle','none')
    hold on;
  end
  view(2)
  grid off
%  lighting phong
  lighting gouraud
  axis equal
  set(gca,'color','w')
  hold off
  