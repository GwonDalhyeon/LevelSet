clc
clear all
close all

l=4/2;
[X,Y]=meshgrid(-l:4/128:l,-l:4/128:l);




% InitialSurface = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\ConsoleApplication1\ConsoleApplication1\initial surface.dat');
% figure
% plot(InitialSurface(:,1),InitialSurface(:,2));
%
figure
initial = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\reinitialzation\reinitialzation\initial level set.dat');
subplot(1,3,1)
contour(initial);
subplot(1,3,2)
surf(X,Y,initial);
subplot(1,3,3)
contour(initial,[0 0])

figure
sign_distance = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\reinitialzation\reinitialzation\sign distance level set.dat');
subplot(1,3,1)
contour(sign_distance);
subplot(1,3,2)
surf(X,Y,sign_distance);
subplot(1,3,3)
contour(sign_distance,[0 0]);


figure
%  hold on
for i=0:1:399
    i
    
    filename = strcat('D:\FE\reinitialization set ',{' '} ,num2str(i),'.dat');
    levelset = importdata(filename{1});
    
    subplot(1,3,1)
    contour(levelset);
    %       hold on
    
    subplot(1,3,2)
    surf(X,Y,levelset);
    %      contour(levelset,[0 0],'color','r')
    subplot(1,3,3)
    contour(levelset,[0 0],'color','r')
    
    
    drawnow
    
    %  clearvars levelset;
    %      hold off
end


figure
%  hold on
for i=0:1:399
    i
    filename = strcat('D:\FE\reinitialization set ',{' '} ,num2str(i),'.dat');
    levelset = importdata(filename{1});
    subplot(1,3,1)
    surf(X,Y,levelset);
    
    filename1 = strcat('D:\FE\Error ',{' '} ,num2str(i),'.dat');
    levelset1 = importdata(filename1{1});
    subplot(1,3,2)
    surf(X,Y,levelset1);
    
    subplot(1,3,3)
    surf(X,Y,sign_distance);
    
    
    drawnow
    
    %  clearvars levelset;
    %      hold off
end



figure
%  hold on
for i=0:1:399
    i
    filename = strcat('D:\FE\reinitialization set ',{' '} ,num2str(i),'.dat');
    levelset = importdata(filename{1});
    subplot(1,3,1)
    surf(X,Y,levelset);
    
    subplot(1,3,2)
    contour(sign_distance,[0 0],'color','r');
    hold on
    contour(levelset,[0 0],'color','b');
    hold off
    
    subplot(1,3,3)
    surf(X,Y,sign_distance);
    
    
    drawnow
    
    %  clearvars levelset;
    %      hold off
end
