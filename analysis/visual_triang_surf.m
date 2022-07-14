clear all
tri_name = '177triangles.txt';
%'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\triangles.txt';
tri = load(tri_name);

vert_name = '177vertices.txt';
%'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\vertices.txt';
vert = load(vert_name);

x = vert(:,1);
y = vert(:,2);
meshgrid(x,y);
z = vert(:,3);
% example:
% [x,y] = meshgrid(1:15,1:15);
% z = peaks(15);
% T = delaunay(x,y);
% trisurf(T,x,y,z)
figure('color','w');
trisurf(tri+1,x,y,z)