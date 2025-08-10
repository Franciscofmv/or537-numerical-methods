% plot examples
function plot_examples
figure(1);
plot_test1;
figure(2);
plot_test2;
figure(3);
plot_test3;
plot_test4;
plot_test5;
plot_test6;
figure(4);
plot_test7;
% figure(5);
% plot_test8;
return

function plot_test1
subplot(2,1,1);
x = -pi:.1:pi;
y1 = 20*sin(x.^2);
y2 = gamma(x);
plot(x,y1,'-ro',x,y2,'-g');
title('function plots');
xlabel('x');
ylabel('y');
legend('20sin(x^2)','gamma');
set(gca,'XTick',-pi:pi/2:pi)
set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
subplot(2,1,2);
x = -2.9:0.1:2.9;
y = randn(10000,1);
hist(y,x);
% h=findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w');
title('histogram');
hold off;
return

function plot_test2
Y = round(rand(5,3)*10);
subplot(2,2,1);
bar(Y,'group');
title('Group');
subplot(2,2,2);
bar(Y,'stack');
title('Stack');
subplot(2,2,3);
barh(Y,'stack');
title('Stack');
subplot(2,2,4);
bar(Y,1.5);
title('Width = 1.5');
hold off;
return

% animation of surface plot
function plot_test3
% divide figure layout into 2x2 quadrants 
%        and place next plot into quardrant 1
subplot(2,2,1); 
Z = peaks; % use peaks function to create 49x49 mesh of Z data
surf(Z); % creat surface plot
axis tight % set axis limits to range of data
set(gca,'nextplot','replacechildren');
for j = 1:20
    surf(sin(2*pi*j/20)*Z,Z); % change plot each time by a sin function
    F(j) = getframe; % return a snapshot of the current axis into F
end
movie(F,2); % Play the movie two times
title('surface animation');
hold off;
return

% animation of line plot
function plot_test4
% another animation
% divide figure layout into 2x2 quadrants and 
%          place next figure into quardrant 2
subplot(2,2,2);
axis tight;
for k = 1:16
    plot(fft(eye(k+16)));
    axis equal;
    M(k) = getframe;
end
movie(M); % play movie once
title('line animation');
hold off;
return

% plot isosurface
function plot_test5
% divide figure layout into 2x2 quadrants 
%        and place next plot into quardrant 1
subplot(2,2,3);
title('3d surface');
[x y z] = meshgrid(1:20,1:20,1:20);
data = sqrt(x.^2 + y.^2 + z.^2);
cdata = smooth3(rand(size(data)),'box',7);
p = patch(isosurface(x,y,z,data,10));
% isonormals(x,y,z,data,p);
isocolors(x,y,z,cdata,p);
set(p,'FaceColor','interp','EdgeColor','none');
view(150,30); 
daspect([1 1 1]); % set equal data aspect ratio for all axes
axis tight;
camlight; % create a light right and up from camera
lighting phong; % set lighting to phong
hold off;
return

% plot isosurface
function plot_test6
subplot(2,2,4);
data = cat(3, [0 .2 0; 0 .3 0; 0 0 0], ...
              [.1 .2 0; 0 1 0; .2 .7 0],...
              [0 .4 .2; .2 .4 0;.1 .1 0]);
data = interp3(data,3,'cubic');
% Draw an isosurface from the volume data and add lights.
% This isosurface uses triangle normals
p1 = patch(isosurface(data,.5),...
'FaceColor','red','EdgeColor','none'); % fill
view(3); 
daspect([1,1,1]); % set equal data aspect ratio for all axes
axis tight; % set axis limits to range of data
camlight(-80,-10); 
lighting phong; 
title('Triangle Normals');
hold off;
return

function plot_test7
axis([0 10 0 10]);
x=[2 4 4 2];
y=[2 2 4 4];
patch(x,y,3);
hold on;
x=[6 8 8 6];
y=[6 6 8 8];
patch(x,y,2);
return

% animation of surface plot and saving in a avi file
function plot_test8
% divide figure layout into 2x2 quadrants 
%        and place next plot into quardrant 1
Z = peaks; % use peaks function to create 49x49 mesh of Z data
surf(Z); % creat surface plot
axis tight % set axis limits to range of data
set(gca,'nextplot','replacechildren');
clear mex; % clear all opened movie files
mov = avifile('plot_example.avi'); % create movie file
mov.Quality = 100;
for j = 1:20
    surf(sin(2*pi*j/20)*Z,Z); % change plot each time by a sin function
    F(j) = getframe; % return a snapshot of the current axis into F
    mov = addframe(mov,F); % add frame to movie file
end
movie(F); % Play the movie just once
title('surface animation');
mov = close(mov); % close movie file
hold off;
return