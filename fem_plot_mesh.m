function [] = fem_plot_mesh(Space);

  figure
  hold
 
  if(strcmpi(Space.e.etype,'triangle'))
    for i = 1:size(Space.t,1)
      x = [Space.x(Space.t(i,1),1:2);
           Space.x(Space.t(i,2),1:2) ];
      plot(x(:,1),x(:,2),'k')

      x = [Space.x(Space.t(i,1),1:2);
           Space.x(Space.t(i,3),1:2) ];
      plot(x(:,1),x(:,2),'k')

      x = [Space.x(Space.t(i,2),1:2);
           Space.x(Space.t(i,3),1:2) ];
      plot(x(:,1),x(:,2),'k')
    end
   
  elseif(strcmpi(Space.e.etype,'quadrilateral'))
    for i = 1:size(Space.t,1)
      x = [Space.x(Space.t(i,1),1:2);
           Space.x(Space.t(i,2),1:2) ];
      plot(x(:,1),x(:,2),'r')

      x = [Space.x(Space.t(i,1),1:2);
           Space.x(Space.t(i,3),1:2) ];
      plot(x(:,1),x(:,2),'r')

      x = [Space.x(Space.t(i,2),1:2);
           Space.x(Space.t(i,4),1:2) ];
      plot(x(:,1),x(:,2),'r')

      x = [Space.x(Space.t(i,3),1:2);
           Space.x(Space.t(i,4),1:2) ];
      plot(x(:,1),x(:,2),'r')
    end
  else
    error('Unknown element type.')
  end
  %plot(Space.x(:,1),Space.x(:,2),'o');

