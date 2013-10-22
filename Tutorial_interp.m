%% ----------------------------------------
%% %% Example file formats
%% ----------------------------------------
%% 1D input format
%{
fid = fopen('interp1D_example.txt','w');
x=0:10;
y = sin(x);

fprintf(fid,'%i %i\n', [size(x,2) -999]);
fprintf(fid, '%f ', x);
fprintf(fid, '\n');
fprintf(fid, '%f ', y);
fprintf(fid, '\n');
fclose(fid);
%% 2D input format
% The first row of Z correspons to y coordinate Y(1) (minimum y)
% The last row of Z corresponds to y coordinate Y(end) maximum y) 
[X,Y] = meshgrid(-3:.25:3);
Z = 1000*peaks(X,Y);
fid = fopen('interp2D_example.txt','w');
fprintf(fid,'%i %i %i\n', [size(X,2) size(Y,1) -999]);
fprintf(fid, '%f ', X(1,:));
fprintf(fid, '\n');
fprintf(fid, '%f ', Y(:,1));
fprintf(fid, '\n');
frmt = '%f';
for i = 2:size(X,2)
    frmt = [frmt ' %f'];
end
frmt = [frmt '\n'];
fprintf(fid, frmt, Z');
fclose(fid);
%% 3D input format
% The first row of V correspons to y coordinate Y(1) (minimum y)
% The last row of V corresponds to y coordinate Y(end) maximum y) 
% The first layer of V corresponds to Z(1) minimum z
% The last layer of V corresponds to Z(end) maximum z

[x,y,z,v] = flow(10);
fid = fopen('interp3D_example.txt','w');
fprintf(fid, '%i %i %i %i\n', [size(x,2) size(y,1), size(z,3) -999]);
fprintf(fid, '%f ', x(1,:,1));
fprintf(fid, '\n');
fprintf(fid, '%f ', y(:,1,1));
fprintf(fid, '\n');
fprintf(fid, '%f ', z(1,1,:));
fprintf(fid, '\n');

frmt = '%f';
for i = 2:size(x,2)
    frmt = [frmt ' %f'];
end
frmt = [frmt '\n'];

for i = 1:size(z,3)
    fprintf(fid, frmt, v(:,:,i)');
end
fclose(fid);

%% 3D input format (varying elevation)
[x,y,z,v] = flow(10);

% change the elevation of z
% The top layer will pass through  the points 
A = [0, -3, 3]; B = [0 3 5]; C=[10 -3 4];
a = (B(2) - A(2))*(C(3) - A(3)) - (C(2) - A(2))*(B(3) - A(3));
b = (B(3) - A(3))*(C(1) - A(1)) - (C(3) - A(3))*(B(1) - A(1));
c = (B(1) - A(1))*(C(2) - A(2)) - (C(1) - A(1))*(B(2) - A(2));
d = -(a*A(1) + b*A(2) + c*A(3));

for i = 1:size(z,3)
    z(:,:,i) = (a*x(:,:,i) + b*y(:,:,i) +d - (i-1)*50)/-c;
end

fid = fopen('interp3D_example1.txt','w');
fprintf(fid, '%i %i %i %i\n', [size(x,2) size(y,1), size(z,3) -999]);
fprintf(fid, '%f ', x(1,:,1));
fprintf(fid, '\n');
fprintf(fid, '%f ', y(:,1,1));
fprintf(fid, '\n');
frmt = '%f';
for i = 2:size(x,2)
    frmt = [frmt ' %f'];
end
frmt = [frmt '\n'];

for i = 1:size(z,3)
    fprintf(fid, frmt, z(z(:,:, size(z,3) - i + 1)');
end
for i = 1:size(z,3)
    fprintf(fid, frmt, v(:,:,i)');
end
fclose(fid);
%}
%% ----------------------------------------
%% %% Examples calling from Matlab/Octave
%% ----------------------------------------
%% 1D
clear; clc;
x=0:10;
y = sin(x);
p = [-1:0.33:11]';
x=x(:);
y=y(:);
v0 = interpND(p,x,[],[],y,-999);

%% 2D
clear; clc;
[X,Y] = meshgrid(-3:.25:3);
Z = 1000*peaks(X,Y);
xx = X(1,:)';
yy = Y(:,1);

[xi yi] = meshgrid(-3.5:0.3:3.5, -3.5:0.3:3.5);
nn = size(xi,1)*size(xi,2);
p = [reshape(xi,nn,1) reshape(yi,nn,1)];
v0 = interpND(p,xx,yy,[],Z,-999);
%% 3D
clear; clc;
load flow_data.mat
xx = x(1,:,1);
yy = y(:,1,1);
zz = reshape(z(1,1,:),length(z(1,1,:)),1);


[xi yi zi] = meshgrid(-0.1:0.4:11, -4.5:0.4:4.5, -4:0.4:4.5);
nn = size(xi,1)*size(xi,2)*size(xi,3);
p = [reshape(xi,nn,1) reshape(yi,nn,1) reshape(zi,nn,1)];
v0 = interpND(p,xx,yy,zz,v,-999);

%% 3D (varying elevation)
clear; clc;
load flow_data.mat
xx = x(1,:,1);
yy = y(:,1,1);
% modify the z elevation
A = [0, -3, 3]; B = [0 3 5]; C=[10 -3 4];
a = (B(2) - A(2))*(C(3) - A(3)) - (C(2) - A(2))*(B(3) - A(3));
b = (B(3) - A(3))*(C(1) - A(1)) - (C(3) - A(3))*(B(1) - A(1));
c = (B(1) - A(1))*(C(2) - A(2)) - (C(1) - A(1))*(B(2) - A(2));
d = -(a*A(1) + b*A(2) + c*A(3));
for i = 1:size(z,3)
    zz(:,:,11-i) = (a*x(:,:,i) + b*y(:,:,i) +d - (i-1)*50)/-c;
end


[xi yi zi] = meshgrid(-1.1:0.4:11, -3.5:0.3:3.5, -4.5:0.4:6.2);
nn = size(xi,1)*size(xi,2)*size(xi,3);
p = [reshape(xi,nn,1) reshape(yi,nn,1) reshape(zi,nn,1)];
v0 = interpND(p,xx,yy,zz,v,-999);



