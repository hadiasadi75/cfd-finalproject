clear 
close all
clc

%% 
Nx = [17 33 65 129]; 
Ny = [17 33 65 129]; 

L = 1;
H = 1;
dx = L./(Nx-1);
dy = H./(Ny-1);

Re = 100; %Reynolds number
nu = 1/Re;
for i=1:4
[u_final,v_final,p_final]=solve_simple(dx(i),dy(i),Nx(i),Ny(i),nu,i);
x = 0:dx(i):L; 
y = 0:dy(i):H; 


%% 
figure(11);hold on
plot(u_final(:,(Ny(i)+1)/2),1-y, 'LineWidth', 1)
figure(111);hold on
plot(y,v_final((Ny(i)+1)/2,:), 'LineWidth', 1)
%% 
x_dom = ((1:Nx(i))-1).*dx(i);
y_dom = 1-((1:Ny(i))-1).*dy(i);
[X,Y] = meshgrid(x_dom,y_dom);
figure;
contourf(X,Y,u_final, 20, 'LineWidth', 1)
colorbar
colormap('jet')
xlabel('x')
ylabel('y')

figure(22);
hold on
quiver(X, Y, u_final, v_final, 5, 'k')
 end

 figure(11); hold on
 plot(u_ghia, y_ghia, 'o', 'LineWidth', 1)
 xlabel('u')
 ylabel('y')
 
 
%     fileID = fopen('lid500.dat','w');
% 
% fprintf(fileID,'%s \n','VARIABLES = "X" , "Y" , "u_final", "v_final", "p_final"');
%  fprintf(fileID,'%s %d %s %d \n',' ZONE I = ',n_points,', J = ',n_points);
%  for i=1:n_points
%     for j=1:n_points 
%         yy=1-(i-1)*h;xx=(j-1)*h;
% fprintf(fileID,'%g %g %g %g %g\r\n',xx,yy, u_final(i,j), v_final(i,j),p_final(i,j));
%     end
%  end
% fclose(fileID);

