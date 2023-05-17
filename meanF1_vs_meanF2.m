clear all;
clc;

filename = '250_4_formula.xls';
File = xlsread(filename);
sheet = 'Sheet1';
xlRange = 'A2:A35035'; %F6
xlRange2 = 'D2:D35035'; % F20

x = xlsread(filename,sheet,xlRange);
y = xlsread(filename,sheet,xlRange2);

prompt1 = 'Enter a floating-point value: ';
a = input(prompt1, 's');
a = str2double(a);

n = size(x,1);
% m = size(rms_vec,2);

temp =0;
for i = 1 : n
    if temp + a > n
        break;
    end

    x_vec(i,:) = x(temp + 1 : temp + a);
    temp = temp +1;
    
end

temp =0;
for i = 1 : n
    if temp + a > n
        break;
    end

    y_vec(i,:) = y(temp + 1 : temp + a);
    temp = temp +1;
    
end

x_mean = mean(x_vec');
y_mean = mean(y_vec');

% We have 2 vectors for mean 
% Code for Video starts from here

filename= 'F6_vs_F20_mean_color.avi';     
vidObj = VideoWriter(filename);
vidObj.FrameRate = 175;
open(vidObj);


h = figure(1);
set(h, 'Position', [50 50 1024 640], 'Color', 'white')
% axis([0 0.4 -100 250])

for i=1:10:length(x_mean)
   
    if i<20496
        plot(x(1:i,1),y(1:i,1),'o','MarkerSize',3,'MarkerFaceColor','b')
    elseif i>=20496 && i<25723
        plot(x(1:20495,1),y(1:20495,1),'o','MarkerSize',3,'MarkerFaceColor','b')
        hold on
        plot(x(20496:i,1),y(20496:i,1),'o','MarkerSize',3,'MarkerFaceColor','r')
        hold off
    else
        plot(x(1:20495,1),y(1:20495,1),'o','MarkerSize',3,'MarkerFaceColor','b')
        hold on
        plot(x(20496:25722,1),y(20496:25722,1),'o','MarkerSize',3,'MarkerFaceColor','r')
        hold on
        plot(x(25723:i,1),y(25723:i,1),'o','MarkerSize',3,'MarkerFaceColor','g')
        hold off
    end
%     hold on
    refreshdata
    drawnow

     axis([0 0.4 -100 250])
    F = getframe(h);
    writeVideo(vidObj,F);
end
% hold off

close(vidObj);






