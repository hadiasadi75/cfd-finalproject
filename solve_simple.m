function [u_final,v_final,p_final]=solve_simple(dx,dy,Nx,Ny,nu,kk)


alpha = 0.01;
alpha_p = 0.8;

u_final(Nx,Ny)=0;
v_final(Nx,Ny)=0;
p_final(Nx,Ny)=1;
u_final(1,:) = 1;


u(Nx+1,Ny)=0;
u_star(Nx+1,Ny)=0;
d_e(Nx+1,Ny)=0;
v(Nx,Ny+1)=0;
v_star(Nx,Ny+1)=0;
d_n(Nx,Ny+1)=0;
p(Nx+1,Ny+1)=1;
p_star(Nx+1,Ny+1)=1;
pp(Nx+1,Ny+1)=0;
b(Nx+1,Ny+1)=0;
u(1,:)=2;

u_new(Nx+1,Ny)=0;
v_new(Nx,Ny+1)=0;
p_new(Nx+1,Ny+1)=1;
u_new(1,:)=2;


error = 1;
iterations = 0;
error_req = 1e-7; 
 figure(kk); 

while error > error_req
   
    for i = 2:Nx
        for j = 2:Ny - 1
            u_E = 0.5*(u(i,j) + u(i,j+1));
            u_W = 0.5*(u(i,j) + u(i,j-1));
            v_N = 0.5*(v(i-1,j) + v(i-1,j+1));
            v_S = 0.5*(v(i,j) + v(i,j+1));
            
            a_E = -0.5*u_E*dy + nu;
            a_W = 0.5*u_W*dy + nu;
            a_N = -0.5*v_N*dx + nu;
            a_S = 0.5*v_S*dx + nu;
            
            a_e = 0.5*u_E*dy - 0.5*u_W*dy + 0.5*v_N*dx - 0.5*v_S*dx + 4*nu;
            
            A_e = -dy;
            d_e(i,j) = A_e/a_e;
            
            u_star_mid = (a_E*u(i,j+1) + a_W*u(i,j-1) + a_N*u(i-1,j) + a_S*u(i+1,j))/a_e + d_e(i,j)*(p(i,j+1) - p(i,j));
            u_star(i,j) = (1-alpha)*u(i,j) + alpha*u_star_mid;
        end
    end
    
    
    u_star(1,:) = 2 - u_star(2,:);
    u_star(Nx + 1,:) = -u_star(Nx,:);
    u_star(2:Nx,1) = 0;
    u_star(2:Nx,Ny) = 0;
    
    
    for i = 2:Nx - 1
        for j = 2:Ny
            u_E = 0.5*(u(i,j) + u(i+1,j));
            u_W = 0.5*(u(i,j-1) + u(i+1,j-1));
            v_N = 0.5*(v(i-1,j) + v(i,j));
            v_S = 0.5*(v(i,j) + v(i+1,j));
            
            a_E = -0.5*u_E*dy + nu;
            a_W = 0.5*u_W*dy + nu;
            a_N = -0.5*v_N*dx + nu;
            a_S = 0.5*v_S*dx + nu;
            
            a_n = 0.5*u_E*dy - 0.5*u_W*dy + 0.5*v_N*dx - 0.5*v_S*dx + 4*nu;
            
            A_n = -dx;
            d_n(i,j) = A_n/a_n;
            
            v_star_mid = (a_E*v(i,j+1) + a_W*v(i,j-1) + a_N*v(i-1,j) + a_S*v(i+1,j))/a_n + d_n(i,j)*(p(i,j) - p(i+1,j));
            v_star(i,j) = (1-alpha)*v(i,j) + alpha*v_star_mid;
        end
    end
    
    % y-momentum eq. - Boundary
    v_star(:,1) = -v_star(:,2);
    v_star(:,Ny + 1) = -v_star(:,Ny);
    v_star(1,2:Ny) = 0;
    v_star(Nx,2:Ny) = 0;
    
    % Zeroing the corrections to begin with
    pp(1:Nx+1,1:Ny+1)=0;
    
    % Continuity equation a.k.a. pressure correction - Interior
    for i = 2:Nx
        for j = 2:Ny
            a_E = -d_e(i,j)*dy;
            a_W = -d_e(i,j-1)*dy;
            a_N = -d_n(i-1,j)*dx;
            a_S = -d_n(i,j)*dx;
            a_P = a_E + a_W + a_N + a_S;
            b(i,j) = -(u_star(i,j) - u_star(i,j-1))*dy + (v_star(i,j) - v_star(i-1,j))*dx;
            
            pp(i,j) = (a_E*pp(i,j+1) + a_W*pp(i,j-1) + a_N*pp(i-1,j) + a_S*pp(i+1,j) + b(i,j))/a_P;
        end
    end
    
    % Correcting the pressure field
    for i = 2:Nx
        for j = 2:Ny
            p_new(i,j) = p(i,j) + alpha_p*pp(i,j);
        end
    end
    
    p_new(1,:) = p_new(2,:);
    p_new(Nx + 1,:) = p_new(Nx,:);
    p_new(:,1) = p_new(:,2);
    p_new(:,Ny + 1) = p_new(:,Ny);
    
    
    for i = 2:Nx
        for j = 2:Ny - 1
            u_new(i,j) = u_star(i,j) + d_e(i,j)*(pp(i,j+1) - pp(i,j));
        end
    end
    
    u_new(1,:) = 2 - u_new(2,:);
    u_new(Nx + 1,:) = -u_new(Nx,:);
    u_new(2:Nx,1) = 0;
    u_new(2:Nx,Ny) = 0;
    
    for i = 2:Nx - 1
        for j = 2:Ny
            v_new(i,j) = v_star(i,j) + d_n(i,j)*(pp(i,j) - pp(i+1,j));
        end
    end
    
    % y-momentum eq. - Boundary
    v_new(:,1) = -v_new(:,2);
    v_new(:,Ny + 1) = -v_new(:,Ny);
    v_new(1,2:Ny) = 0;
    v_new(Nx,2:Ny) = 0;
            
    
    % Continuity residual as error measure
    error = 0;
    for i = 2:Nx
        for j = 2:Ny
            error = error + abs(b(i,j));
        end
    end
    
   
    if(rem(iterations, 1000)) == 0
       figure(kk);
       semilogy(iterations, error, '-ko')
       hold on
       xlabel('Iterations')
       ylabel('Residual Error')
    end
    
    % Finishing the iteration
    u = u_new;
    v = v_new;
    p = p_new;
    iterations = iterations + 1;
    
end


for i = 1:Nx
    for j = 1:Nx
        u_final(i,j) = 0.5*(u(i,j) + u(i+1,j));
        v_final(i,j) = 0.5*(v(i,j) + v(i,j+1));
        p_final(i,j) = 0.25*(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1));
    end
end


    

