%%
[data,fs] = audioread("1-s2.0-S2666386421000722-mmc4.mp3");
% data- input data, fs-sampling data
% size(data)

N = length(data);
t = (0:N-1)/fs;



%%
window_size = '250';
window_overlap = '50';
prompt1 = 'Enter a floating-point value: ';
x = input(prompt1, 's');
x = str2double(x);
filename = '250_4_formula.xls';
sheet = 1;
%%
figure(250);

subplot(4,1,1);
plot(t,data');

title('AE Intensity vs time');
xlabel('t (s)');
ylabel('A (a.u.)');
ax = gca;
annotation('textbox',[ax.Position(1)-0.1,ax.Position(2)+ax.Position(4)-0.02,0.09,0.09],'String',['Window size: ' num2str(window_size) 'ms'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');
annotation('textbox',[ax.Position(1) + 0.6, ax.Position(2) + ax.Position(4) - 0.02, 0.09, 0.09], 'String', ['Window overlap: ' num2str(window_overlap) '%'], 'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');

hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

xlRange = 'A2:A35035';

rms_vec = xlsread(filename,sheet,xlRange);
n = size(rms_vec,1);
m = size(rms_vec,2);

tn = linspace(0,t(end),n);

subplot(4,1,2);
plot(tn,rms_vec');
xlabel('t (s)');
ylabel('F6');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;
temp =0;
for i = 1 : n
    if temp + x > n
        break;
    end

    rms(i,:) = rms_vec(temp + 1 : temp + x);
    temp = temp +1;
    
end  

r = size(rms,1);
t_x = linspace(0,t(end),r);
rms_std = std(rms');
rms_mean = mean(rms');

subplot(4,1,3);
plot(t_x,rms_std);
xlabel('t (s)');
ylabel('Std of F6');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

subplot(4,1,4);
plot(t_x,rms_mean);
xlabel('t (s)');
ylabel('Mean of F6');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;
%%
xlRange = 'B2:B35035';

f53_vec = xlsread(filename,sheet,xlRange);
n = size(f53_vec,1);
m = size(f53_vec,2);

tn = linspace(0,t(end),n);

figure(251);
subplot(4,1,1);
plot(t,data');
title('AE Intensity vs time');
xlabel('t (s)');
ylabel('A (a.u.)');
ax = gca;
annotation('textbox',[ax.Position(1)-0.1,ax.Position(2)+ax.Position(4)-0.02,0.09,0.09],'String',['Window size: ' num2str(window_size) 'ms'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');
annotation('textbox',[ax.Position(1) + 0.6, ax.Position(2) + ax.Position(4) - 0.02, 0.09, 0.09], 'String', ['Window overlap: ' num2str(window_overlap) '%'], 'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');

hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;



subplot(4,1,2);
plot(tn,f53_vec');
xlabel('t (s)');
ylabel('F53');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

temp =0;
for i = 1 : n
    if temp + x > n
        break;
    end

    f53(i,:) = f53_vec(temp + 1 : temp + x);
    temp = temp +1;
    
end  

r = size(f53,1);
t_x = linspace(0,t(end),r);
f53_std = std(f53');
f53_mean = mean(f53');

subplot(4,1,3);
plot(t_x,f53_std);
xlabel('t (s)');
ylabel('Std of F53');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

subplot(4,1,4);
plot(t_x,f53_mean);
xlabel('t (s)');
ylabel('Mean of F53');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

%%
xlRange = 'C2:C35035';

f59_vec = xlsread(filename,sheet,xlRange);
n = size(f59_vec,1);
m = size(f59_vec,2);

tn = linspace(0,t(end),n);

figure(252);
subplot(4,1,1);
plot(t,data');

title('AE Intensity vs time');
xlabel('t (s)');
ylabel('A (a.u.)');
ax = gca;
annotation('textbox',[ax.Position(1)-0.1,ax.Position(2)+ax.Position(4)-0.02,0.09,0.09],'String',['Window size: ' num2str(window_size) 'ms'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');
annotation('textbox',[ax.Position(1) + 0.6, ax.Position(2) + ax.Position(4) - 0.02, 0.09, 0.09], 'String', ['Window overlap: ' num2str(window_overlap) '%'], 'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');

hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;


subplot(4,1,2);
plot(tn,f59_vec');
xlabel('t (s)');
ylabel('F59');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

temp =0;
for i = 1 : n
    if temp + x > n
        break;
    end

    f59(i,:) = f59_vec(temp + 1 : temp + x);
    temp = temp +1;
    
end  

r = size(f59,1);
t_x = linspace(0,t(end),r);
f59_std = std(f59');
f59_mean = mean(f59');

subplot(4,1,3);
plot(t_x,f59_std);
xlabel('t (s)');
ylabel('Std of F59');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

subplot(4,1,4);
plot(t_x,f59_mean);
xlabel('t (s)');
ylabel('Mean of F59');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

%%
xlRange = 'D2:D35035';

f20_vec = xlsread(filename,sheet,xlRange);
n = size(f20_vec,1);
m = size(f20_vec,2);

tn = linspace(0,t(end),n);

figure(253);
subplot(4,1,1);
plot(t,data');
title('AE Intensity vs time');
xlabel('t (s)');
ylabel('A (a.u.)');
ax = gca;
annotation('textbox',[ax.Position(1)-0.1,ax.Position(2)+ax.Position(4)-0.02,0.09,0.09],'String',['Window size: ' num2str(window_size) 'ms'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');
annotation('textbox',[ax.Position(1) + 0.6, ax.Position(2) + ax.Position(4) - 0.02, 0.09, 0.09], 'String', ['Window overlap: ' num2str(window_overlap) '%'], 'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');

hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

subplot(4,1,2);
plot(tn,f20_vec');
xlabel('t (s)');
ylabel('F20');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

temp =0;
for i = 1 : n
    if temp + x > n
        break;
    end

    f20(i,:) = f20_vec(temp + 1 : temp + x);
    temp = temp +1;
    
end  

r = size(f20,1);
t_x = linspace(0,t(end),r);
f20_std = std(f20');
f20_mean = mean(f20');

subplot(4,1,3);
plot(t_x,f20_std);
xlabel('t (s)');
ylabel('Std of F20');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

subplot(4,1,4);
plot(t_x,f20_mean);
xlabel('t (s)');
ylabel('Mean of F20');
hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;



