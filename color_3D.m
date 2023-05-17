clear all;
clc;

filename = '250_4_formula.xls';
sheet = 1;

xlRange = 'A2:A35035';
ylRange = 'B2:B35035';
zlRange = 'C2:C35035';
x = xlsread(filename,sheet,xlRange);
y = xlsread(filename,sheet,ylRange);
z = xlsread(filename,sheet,zlRange);



filename= 'F6_F53_F59.avi';     
vidObj = VideoWriter(filename);
vidObj.FrameRate = 175;
open(vidObj);


h = figure(1);
set(h, 'Position', [50 50 1024 640], 'Color', 'white')
% axis([0 0.4 -100 250])

for i=1:10:length(x)
   
    if i<20496
        plot3(x(1:i,1),y(1:i,1),z(1:i,1),'o','MarkerSize',3,'MarkerFaceColor','b')
        grid on % turn on the grid
    elseif i>=20496 && i<25723
        plot3(x(1:20495,1),y(1:20495,1),z(1:20495,1),'o','MarkerSize',3,'MarkerFaceColor','b')
        grid on % turn on the grid
        hold on
        plot3(x(20496:i,1),y(20496:i,1),z(20496:i,1),'o','MarkerSize',3,'MarkerFaceColor','r')
        grid on % turn on the grid
        hold off
    else
        plot3(x(1:20495,1),y(1:20495,1),z(1:20495,1),'o','MarkerSize',3,'MarkerFaceColor','b')
        grid on % turn on the grid
        hold on
        plot3(x(20496:25722,1),y(20496:25722,1),z(20496:25722,1),'o','MarkerSize',3,'MarkerFaceColor','r')
        grid on % turn on the grid
        hold on
        plot3(x(25723:i,1),y(25723:i,1),z(25723:i,1),'o','MarkerSize',3,'MarkerFaceColor','g')
        grid on % turn on the grid
        hold off
    end
%     hold on
    refreshdata
    drawnow
    
    xlabel('F6')
    ylabel('F53')
    zlabel('F59')
    axis([0 0.4 -1 0 -0.08 0.1])
    F = getframe(h);
    writeVideo(vidObj,F);
end
% hold off

close(vidObj);