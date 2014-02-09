clc
clear all
close all

l=4/2;
[X,Y]=meshgrid(-l:4/150:l,-l:4/150:l);




% InitialSurface = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\ConsoleApplication1\ConsoleApplication1\initial surface.dat');
% figure
% plot(InitialSurface(:,1),InitialSurface(:,2));
% 
levelset = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\reinitialzation\reinitialzation\initial level set.dat');
figure
 surf(X,Y,levelset);
 

 figure
 contour(levelset,[0 0]);
 
 figure
 contour(levelset);
 
 figure
%  hold on
 for i=0:1:1999
     i
    
     filename = strcat('D:\FE\reinitialization set ',{' '} ,num2str(i),'.dat');
     levelset = importdata(filename{1});
     
%      subplot(1,2,1)
%      contour(levelset);
%       hold on

%      subplot(1,2,2)
     surf(X,Y,levelset);
%      contour(levelset,[0 0],'color','r')    
%      
     
     drawnow
     
     %  clearvars levelset;
%      hold off
 end
