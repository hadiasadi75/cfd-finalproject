clear all
close all
clc

%% 

Nx = 51 ; 
L = 1; 
h = L/(Nx - 1); 
x = 0:h:L; 
y = 0:h:L; 
z = 0:h:L; 
dt = 0.001; 
Re = 100; % Reynolds number
delta = 1.5; % Artificial compressibility


%% Initializing the variables

% Final collocated variables

u_final(Nx,Nx,Nx) = 0;
v_final(Nx,Nx,Nx) = 0;
w_final(Nx,Nx,Nx) = 0;
p_final(Nx,Nx,Nx) = 0;
u_final(1,:,:) = 1;

% Staggered variables
u(Nx+1, Nx, Nx+1) = 0;
v(Nx, Nx+1, Nx+1) = 0;
w(Nx+1, Nx+1, Nx) = 0;
p(Nx+1 , Nx+1, Nx+1) = 1;
u(1,:,:) = 2;


u_new(Nx+1, Nx, Nx+1) = 0;
v_new(Nx, Nx+1, Nx+1) = 0;
w_new(Nx+1, Nx+1, Nx) = 0;
p_new(Nx+1 , Nx+1, Nx+1) = 1;
u_new(1,:,:) = 2;


error = 1;
iterations = 0;
error_req = 1e-5; 
figure(1); 

while error > error_req
   
    for k=2:Nx
    for i = 2:Nx
        for j = 2:Nx - 1
            pressure = -(p(i,j+1,k) - p(i,j,k))/h;
            diffusion = (1/Re)*((u(i+1,j,k) - 2*u(i,j,k) + u(i-1,j,k))/(h*h) +...
                (u(i,j+1,k) - 2*u(i,j,k) + u(i,j-1,k))/(h*h)+...
                (u(i,j,k+1) - 2*u(i,j,k) + u(i,j,k-1))/(h*h) );
            advection_x = ((0.5*(u(i,j,k)+u(i,j+1,k)))^2 - (0.5*(u(i,j,k)+u(i,j-1,k)))^2)/h;
            advection_y = ((0.25*(u(i,j,k)+u(i-1,j,k))*(v(i-1,j,k)+v(i-1,j+1,k)))...
                         - (0.25*(u(i,j,k)+u(i+1,j,k))*(v(i,j,k)+v(i,j+1,k))))/h;
            advection_z = ((0.25*(u(i,j,k)+u(i,j,k+1))*(w(i,j,k)+w(i,j+1,k)))...
                         - (0.25*(u(i,j,k)+u(i,j,k-1))*(w(i,j,k-1)+w(i,j+1,k-1))))/h;        
            u_new(i,j,k) = u(i,j,k) + dt*(diffusion - advection_x - advection_y- advection_z+ pressure);
        end
    end
    end
   
    u_new(:,1,:)=0;
    u_new(:,Nx,:)=0;
    u_new(:,:,1)=-u_new(:,:,2);
    u_new(:,:,Nx)=-u_new(:,:,Nx+1);
    u_new(Nx,:,:)=-u_new(Nx+1,:,:);
     u_new(1,:,:)=2-u_new(2,:,:);
    

    for k=2:Nx
    for i = 2:Nx - 1
        for j = 2:Nx
            pressure = -(p(i,j,k) - p(i+1,j,k))/h;
            diffusion = (1/Re)*((v(i+1,j,k) - 2*v(i,j,k) + v(i-1,j,k))/(h*h) +...
                (v(i,j+1,k) - 2*v(i,j,k) + v(i,j-1,k))/(h*h)+...
                (v(i,j,k+1) - 2*v(i,j,k) + v(i,j,k-1))/(h*h) );
            advection_y = ((0.5*(v(i,j,k)+v(i-1,j,k)))^2 - (0.5*(v(i,j,k)+v(i+1,j,k)))^2)/h;
            advection_x = ((0.25*(u(i,j,k)+u(i-1,j,k))*(v(i,j,k)+v(i,j+1,k))) - (0.25*(u(i,j-1,k)+u(i-1,j-1,k))*(v(i,j,k)+v(i,j-1,k))))/h;
            advection_z = ((0.25*(w(i,j,k)+w(i-1,j,k))*(v(i,j,k)+v(i,j,k+1))) - (0.25*(w(i,j,k-1)+w(i-1,j,k-1))*(v(i,j,k)+v(i,j,k-1))))/h;
            v_new(i,j,k) = v(i,j,k) + dt*(diffusion - advection_x - advection_y-advection_z + pressure);
        end
    end
    end
    
     v_new(:,1,:) = -v_new(:,2,:);
    v_new(:,Nx + 1,:) = -v_new(:,Nx,:);
    v_new(1,:,:) = 0;
    v_new(Nx,:,:) = 0;
     v_new(:,:,1)=- v_new(:,:,2);
    v_new(:,:,Nx + 1) = -v_new(:,:,Nx );
    
    for k=2:Nx-1
    for i = 2:Nx 
        for j = 2:Nx
            pressure = -(p(i,j,k+1) - p(i,j,k))/h;
            diffusion = (1/Re)*((w(i+1,j,k) - 2*w(i,j,k) + w(i-1,j,k))/(h*h) +...
                (w(i,j+1,k) - 2*w(i,j,k) + w(i,j-1,k))/(h*h)+...
                (w(i,j,k+1) - 2*w(i,j,k) + w(i,j,k-1))/(h*h) );
            advection_z = ((0.5*(w(i,j,k)+w(i,j,k-1)))^2 - (0.5*(w(i,j,k)+w(i,j,k+1)))^2)/h;
            advection_x = ((0.25*(u(i,j,k)+u(i,j,k+1))*(w(i,j,k)+w(i,j+1,k))) - (0.25*(u(i,j-1,k)+u(i,j-1,k+1))*(w(i,j,k)+w(i,j-1,k))))/h;
            advection_y = ((0.25*(w(i,j,k)+w(i-1,j,k))*(v(i-1,j,k)+v(i-1,j,k+1))) - (0.25*(w(i,j,k)+w(i+1,j,k))*(v(i,j,k)+v(i,j,k+1))))/h;
            w_new(i,j,k) = w(i,j,k) + dt*(diffusion - advection_x - advection_y-advection_z + pressure);
        end
    end
    end
    
   
     w_new(:,:,1) = 0;
    w_new(:,:,Nx + 1) = 0;
    w_new(1,:,:) = -w_new(2,:,:);
    w_new(Nx,:,:) = -w_new(Nx+1,:,:);
     w_new(:,1,:)=- w_new(:,2,:);
    w_new(:,Nx + 1,:) = -w_new(:,Nx,: );
    
    
    for k = 2:Nx
    for i = 2:Nx
        for j = 2:Nx
            p_new(i,j,k) = p(i,j,k) - delta*dt*((u(i,j,k) - u(i,j-1,k) + v(i-1,j,k) - v(i,j,k) + w(i,j,k) - w(i,j,k-1)))/h;
        end
    end
    end
    
    
    
    p_new(1,:,:) = p_new(2,:,:);
    p_new(Nx + 1,:,:) = p_new(Nx,:,:);
    p_new(:,1,:) = p_new(:,2,:);
    p_new(:,Nx + 1,:) = p_new(:,Nx,:);
    p_new(:,:,1) = p_new(:,:,2);
     p_new(:,:,Nx + 1) = p_new(:,:,Nx);
    
     error = 0;
  for k = 2:Nx - 1
    for i = 2:Nx - 1
        for j = 2:Nx - 1
            error = error + abs((u_new(i,j,k) - u_new(i,j-1,k) + v_new(i-1,j,k) - v_new(i,j,k) + w_new(i,j,k) - w_new(i,j,k-1))/h);
        end
    end
  end
    
    if(rem(iterations, 5000)) == 0
       figure(1);
       semilogy(iterations, error, '-ko')
       hold on
       xlabel('Iterations')
       ylabel('Residual Error')
    end
    
   
    u = u_new;
    v = v_new;
      w = w_new;
    p = p_new;
    iterations = iterations + 1;
end


for k = 1:Nx
for i = 1:Nx
    for j = 1:Nx
        u_final(i,j,k) = 0.25*(u(i,j,k) + u(i+1,j,k)+u(i,j,k+1) + u(i+1,j,k+1));
        v_final(i,j,k) = 0.25*(v(i,j,k) + v(i,j+1,k)+v(i,j,k+1) + v(i,j+1,k+1));
        w_final(i,j,k) = 0.25*(w(i+1,j+1,k) + w(i+1,j,k)+w(i,j,k) + w(i,j+1,k));
       
        p_final(i,j,k) = 0.125*(p(i,j,k) + p(i+1,j,k)+p(i,j+1,k) + p(i+1,j+1,k) +...
                              p(i,j,k+1) + p(i+1,j,k+1)+p(i,j+1,k+1) + p(i+1,j+1,k+1));
    end
end
end
  fileID = fopen('lid3d.dat','w');

fprintf(fileID,'%s \n','VARIABLES = "X" , "Y" , "z", "u_final", "v_final" ,"w_final","p_final"');
 fprintf(fileID,'%s %d %s %d %s %d\n',' ZONE I = ',Nx,', J = ',Nx,', k = ',Nx);
 for i=1:Nx
    for j=1:Nx 
        for k=1:Nx 
        yy=1-(i-1)*h;xx=(j-1)*h;zz=(k-1)*h;
fprintf(fileID,'%g %g %g %g %g %g %g\r\n',xx,yy,zz, u_final(i,j,k), v_final(i,j,k), w_final(i,j,k), p_final(i,j,k));
        end
    end
 end
fclose(fileID);

