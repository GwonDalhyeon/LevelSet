clc
clear all
close all

l=4/2;
[X,Y]=meshgrid(-l:4/150:l,-l:4/150:l);




% InitialSurface = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\ConsoleApplication1\ConsoleApplication1\initial surface.dat');
% figure
% plot(InitialSurface(:,1),InitialSurface(:,2));
%
figure('units','normalized','outerposition',[0 0 1 1])
initial = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\reinitialzation\reinitialzation\initial level set.dat');
subplot(1,3,1)
contour(initial);
subplot(1,3,2)
surf(X,Y,initial);
subplot(1,3,3)
contour(initial,[0 0])

figure('units','normalized','outerposition',[0 0 1 1])
sign_distance = importdata('C:\Users\Gwon Dalhyeon\Documents\Visual Studio 2012\Projects\Level Set\reinitialzation\reinitialzation\sign distance function.dat');
subplot(1,3,1)
contour(sign_distance);
subplot(1,3,2)
surf(X,Y,sign_distance);
subplot(1,3,3)
contour(sign_distance,[0 0]);


figure('units','normalized','outerposition',[0 0 1 1])
%  hold on
for i=0:1:399
    i
    %     TVD RK
    filename = strcat('D:\TVD RK\reinitialization set ',{' '} ,num2str(i),'.dat');
    levelset = importdata(filename{1});
    subplot(2,3,1)
    figure
    surf(X,Y,levelset);
%     title(strcat('time',num2str(i)));
    
    filename1 = strcat('D:\TVD RK\Error ',{' '} ,num2str(i),'.dat');
    filename2 = strcat('D:\TVD RK\L0 ',{' '} ,num2str(i),'.dat');
    filename3 = strcat('D:\TVD RK\L1 ',{' '} ,num2str(i),'.dat');
    levelset1 = importdata(filename1{1});
    l0=importdata(filename2{1});
    l1=importdata(filename2{1});
    subplot(2,3,4)
    surf(X,Y,levelset1);
%     title(strcat('L0 error =', num2str(l0), ', L0 error =', num2str(l1)) );
    
    %     FE
    filename = strcat('D:\FE\reinitialization set ',{' '} ,num2str(i),'.dat');
    levelset = importdata(filename{1});
    subplot(2,3,2)
    surf(X,Y,levelset);
    
    filename1 = strcat('D:\FE\Error ',{' '} ,num2str(i),'.dat');
    filename2 = strcat('D:\FE\L0 ',{' '} ,num2str(i),'.dat');
    filename3 = strcat('D:\FE\L1 ',{' '} ,num2str(i),'.dat');
    levelset1 = importdata(filename1{1});
    l0=importdata(filename2{1});
    l1=importdata(filename2{1});
    subplot(2,3,5)
    surf(X,Y,levelset1);
%     title(strcat('L0 error =', num2str(l0), ', L0 error =', num2str(l1)) );
    
    %     GS
    filename = strcat('D:\Gauss Seidel\reinitialization set ',{' '} ,num2str(i),'.dat');
    levelset = importdata(filename{1});
    subplot(2,3,3)
    surf(X,Y,levelset);
    
    filename1 = strcat('D:\Gauss Seidel\Error ',{' '} ,num2str(i),'.dat');
    filename2 = strcat('D:\Gauss Seidel\L0 ',{' '} ,num2str(i),'.dat');
    filename3 = strcat('D:\Gauss Seidel\L1 ',{' '} ,num2str(i),'.dat');
    levelset1 = importdata(filename1{1});
    l0=importdata(filename2{1});
    l1=importdata(filename2{1});
    subplot(2,3,6)
    surf(X,Y,levelset1);
%     title(strcat('L0 error =', num2str(l0), ', L0 error =', num2str(l1)) );
   
    drawnow
    
    %  clearvars levelset;
    %      hold off
end