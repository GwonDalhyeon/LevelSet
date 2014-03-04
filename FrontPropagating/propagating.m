clc
clear all
close all


[X,Y]=meshgrid(-0.5:1/600:0.5,-0.5:1/600:0.5);

InitialSurface = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\ConsoleApplication1\ConsoleApplication1\initial surface.dat');
levelset = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\ConsoleApplication1\ConsoleApplication1\initial level set.dat');


figure
plot(InitialSurface(:,1),InitialSurface(:,2));


figure
 surf(X,Y,levelset);

 figure
 contour(levelset,[0 0]);
 
 figure
 contour(levelset);
 
 
%  levelset10 = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\ConsoleApplication1\ConsoleApplication1\level set 10.dat');
%  figure
%  contour(levelset10);
%  
%  levelset20 = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\ConsoleApplication1\ConsoleApplication1\level set 20.dat');
%  figure
% contour(levelset20);
 

 levelset700 = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\ConsoleApplication1\ConsoleApplication1\level set 700.dat');
 figure
 contour(levelset700);
 
 figure
 for i=0:100:3999
     i
 filename = strcat('D:\level set',{' '} ,num2str(i),'.dat');
 levelset = importdata(filename{1});
 contour(levelset,[0 0])
 drawnow
 clearvars levelset;
 end
 