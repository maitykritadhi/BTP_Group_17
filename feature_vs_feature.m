clear all;
clc;

load('prevData.mat')

filename = '250_4_formula.xls';
File = xlsread(filename);
sheet = 'Sheet1';
xlRange = 'A2:A35035'; %F6
xlRange2 = 'D2:D35035'; % F20


x = xlsread(filename,sheet,xlRange);
y = xlsread(filename,sheet,xlRange2);

filename= 'F6_vs_F20_f2.avi';     
vidObj = VideoWriter(filename);
vidObj.FrameRate = 1750;
open(vidObj);

for i=1:35034
    h = figure(1);
    set(h, 'Position', [50 50 1024 640], 'Color', 'white')
%     set(h, 'Position', [50 50 1024 640], 'Color', 'white')
%     subplot(8,2,[1 10]);
    plot(x(1:i,1),y(1:i,1),'o','MarkerSize',3,'MarkerFaceColor','b')
%     hold on;
     axis([0 0.4 -100 250])
% pause(0.1)
    F = getframe(h);
    writeVideo(vidObj,F);
end

close all;
clear all;
close(vidObj);
