%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                %
%  main solver script for a finite element analysis on the VALVE FLOW PROJECT,   % 
%  AFEM final coursework.                                                        % 
%  by Cesare Magnetti                                                            %
%  Mar 2020                                                                      %
%                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  close all; clear; clc;
  
  %max order term in weakform is v*grad(v)*w, recalling that v and w are
  %approximated by quadratic interpolators: v= 2nd order, grad(v) = 1st
  %order, w = 2nd order --> 5th order polynomial is perfectly integrated
  %using 5 gauss points.
  quad=5; 
  etype = 'triangle';

  nonlinear_tol = 1e-6; % setting the nonlinear tolerance
  residual_vector_norm = 1; % initializing the residual vector norm

% Making Domain function space ...
  load('jetflow3.mat');
  Omega.e=fem_get_basis(Omega.p, quad, etype);

% Making vel variable function space ...
  vel = Omega2;
  vel.name = 'vel'; % field name ...
  vel.e=fem_get_basis(vel.p, quad, etype); 

% Setting our initial guess for velocity to be zero ...
  vel.u = zeros(vel.dm * size(vel.x,1), 1);
  n_vel = size(vel.u,1); 

% Making variable function space ...
  pres = Omega;
  pres.name = 'pres'; % field name ...
  pres.dm = 1; % scalar field ...
  pres.e=fem_get_basis(pres.p, quad, etype); 

  % Setting our initial guess for pressre to be zero ...
  pres.u = zeros(size(pres.x,1),1); 
  n_pres = size(pres.u,1);

% Starting the Newton-Raphson iteration ...
disp(['        A second for loop was added to dampen the effects of the boundary conditions,']);
disp(['        this should stabilize the convergence of Newton Rapson method.']);
tic;
for k = 0.05:0.05:1
    iter = 0; % initializing the interation counter
    while (iter == 0) || (residual_vector_norm > nonlinear_tol) 
      % Making a FEM object storing unknown variables (must do this every iteration) ...
        vars.vel = vel;
        vars.pres = pres;

      % Making the residual vector and applying boundary conditions (velocity) ...
        R1 = fem_assemble_block_residual(@navier_stokes_eqn, Omega, vel, vars);
        
      % variable vel.b contains all the necessary information to apply B.C.
      % the first and last column respectevely report the element id and
      % the type of boundary, in this case we have 4 boundaries numberd 1
      % to 4. The number inner columns vary depending on the order (p)
      % chosen, and it reports the singular nodes of the element that are
      % on the boundary. Note that some elements, and even some nodes, can
      % lie on more than 1 boundary (i.e. at boundary interfaces).
        for i = 1:size(vel.b,1)
          switch vel.b(i,end)
              case 1
                  for ii = 2:size(vel.b,2)-1
                    nn = vel.dm * (vel.b(i,ii) - 1);
                    
                    %contant inflow
                    R1(nn+1:nn+2) = -[-k*200 - vel.u(nn+1); 0 - vel.u(nn+2)];
                  end
              case 2
                  %do nothing, neumann condition disappears in weak form so
                  %we can skip elements on this boundary
                  continue;
              case 3
                  %dirichlet depends on value of x2
                  for ii = 2:size(vel.b,2)-1
                      if vel.x(vel.b(i,ii),2)>0
                          nn = vel.dm * (vel.b(i,ii) - 1);
                          R1(nn+1:nn+2) = -[ 0 - vel.u(nn+1); k*10 - vel.u(nn+2)];
                      else
                          nn = vel.dm * (vel.b(i,ii) - 1);
                          R1(nn+1:nn+2) = -[ 0 - vel.u(nn+1); -k*10 - vel.u(nn+2)];
                      end
                  end   
              case 4
                  for ii = 2:size(vel.b,2)-1
                    nn = vel.dm * (vel.b(i,ii) - 1);
                    R1(nn+1:nn+2) = -[0 - vel.u(nn+1); 0 - vel.u(nn+2)];
                  end
          end

        end

      % Making the residual vector and applying boundary conditions (on species B) ...
        R2 = fem_assemble_block_residual(@darcy_mass_eqn, Omega, pres, vars);

      % Making the global residual vector + computing the norm ...
        R = [R1; R2];
        residual_vector_norm = norm(R,2);
        disp(['Current residual (iter=' num2str(iter) ' , k=' num2str(k) '): ' num2str(residual_vector_norm)])
      % escape condition
        if(residual_vector_norm < nonlinear_tol)
          continue
        end

      % Creating the discrete matrix operators corresponding to operators a, b, and c
        disp(['        constructing the Jacobian blocks ...'])
        A  = fem_assemble_block_matrix_perturbation(@navier_stokes_eqn, Omega,vel, vel , vars); 
        B  = fem_assemble_block_matrix_perturbation(@navier_stokes_eqn, Omega, vel, pres , vars); 
        C  = fem_assemble_block_matrix_perturbation(@darcy_mass_eqn, Omega, pres , vel, vars); 
        D = sparse(size(pres.u,1),size(pres.u,1)); %R2 independent of pressure, derivatives = 0


      % Editing block matrices for dirichlet conditions (we dont want to modify solution at boundaries) ...
        for i = 1:size(vel.b,1)
          if(vel.b(i,end) == 2)
              continue;
          end

          for ii = 2:size(vel.b,2)-1
              nn = vel.dm * (vel.b(i,ii) - 1);
              A(nn+1:nn+2,:) = 0;
              A(nn+1:nn+2, nn+1:nn+2) = eye(2);
              B(nn+1:nn+2,:) = 0;
          end
        end

      % Composing the Jacobian from our Jacobian blocks ...
        disp(['        assembly of the Global Jacobian ...'])
        J = [A B; C D];

      % Apply Newton Raphson ...
        disp(['        solving for the NR update ...'])
        U = J \ R;
        vel.u(1:n_vel) = vel.u(1:n_vel) - U(1:n_vel);
        pres.u(1:n_pres) = pres.u(1:n_pres) - U(n_vel+1:end);
        
      % Update the iteration counter ...
        iter = iter + 1;
        
        disp(['  ']) % skip a line so output is prettier
    end
   %save in between results if needed (make sure that folder jetflowX
   %exists before)
   %save(sprintf('jetflow1/results%0.2f.mat',k), 'pres', 'vel'); 
end
% evaluate convergence time
elapsed_time = toc;
%save results
save('results3.mat', 'pres', 'vel');
%print elapsed time
disp(['convergence time: ' num2str(floor((elapsed_time)/60)) ' min    ' round(num2str(mod(elapsed_time,60))) ' sec']);