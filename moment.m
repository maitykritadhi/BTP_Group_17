%% Formation of basic vectors
B = sortrows(b,1);
t = B(:,1);
data = B(:,2);
heat_flux = B(:,3);
N = length(data);

[maxVal, maxInd] = max(heat_flux); % find the maximum value and its index in the third column of B
t_CHF = t(maxInd);

tem = 0;
for i=1 : N
    tem = tem +1;
    if(t(i)>3000)
        break;
    end
end

% plot(B(:,1), B(:,2));
% hold on;
% plot(t_CHF*ones(size(t)),data', '-r') %CHF Line indicator
% hold off;

% time= linspace(t(1), t(end), N);
% plot(time,t);
%% making of vector frames 
t_max = max(t);
fs = round(N/t_max);
f_d = 2.5;

% Convert the window size to samples
f_size = f_d * fs;

% Calculate the number of frames
n = length(data);

% Extract the frames
temp = 0;
for i = 1 : tem
    frames(i,:) = data(temp + 1 : temp + f_size);
    temp = temp + f_size;
    if(temp + f_size > tem)
        break;
    end
end 
% writematrix(frames,'Frames.xls','Sheet', 1)
[num_rows, num_columns] = size(frames);
a = (3025.27/3443.83);
c = floor(3000/a);
tn = linspace(0,c,num_rows);
tn = tn.*a;
%% Calculating moment vector

% F20 
mean_frames = mean(frames');
% plot(tn, mean_frames);
% hold on;
% plot(t_CHF*ones(size(tn)),mean_frames', '-r') %CHF Line indicator
% hold off;

for i = 1:num_rows
    frames1(i,:) = frames(i,:)-mean_frames(i);
end
frames1 = frames1.^2;
[rows,col] = size(frames1);

moment_5 = sum(frames1');

ten = 1.0e+14;

% Declaring xl range +

% Define the range variable
myrange = 'K5:K16';
%%
% Threshold 0.05
cnt_1 = 0;
thr_1 = 0.05*ten; % threshold
for i=1:size(moment_5')
    cnt_1 = cnt_1 + 1;
    if(moment_5(i)>= thr_1)
        break;
    end
end
t_thr_1 = tn(cnt_1);

[d_1, i_1] = min(abs(t - t_thr_1));
q_thr_1 = heat_flux(i_1);

pct_1 = 100*q_thr_1/max(heat_flux);

% Threshold 0.1
cnt_2 = 0;
thr_2 = 0.1*ten; % threshold
for i=1:size(moment_5')
    cnt_2 = cnt_2 + 1;
    if(moment_5(i)>= thr_2)
        break;
    end
end
t_thr_2 = tn(cnt_2);

[d_2, i_2] = min(abs(t - t_thr_2));
q_thr_2 = heat_flux(i_2);

pct_2 = 100*q_thr_2/max(heat_flux);

% Threshold 0.5
cnt_3 = 0;
thr_3 = 0.5*ten; % threshold
for i=1:size(moment_5')
    cnt_3 = cnt_3 + 1;
    if(moment_5(i)>= thr_3)
        break;
    end
end
t_thr_3 = tn(cnt_3);

[d_3, i_3] = min(abs(t - t_thr_3));
q_thr_3 = heat_flux(i_3);

pct_3 = 100*q_thr_3/max(heat_flux);

% Threshold 1
cnt_4 = 0;
thr_4 = 1*ten; % threshold
for i=1:size(moment_5')
    cnt_4 = cnt_4 + 1;
    if(moment_5(i)>= thr_4)
        break;
    end
end
t_thr_4 = tn(cnt_4);

[d_4, i_4] = min(abs(t - t_thr_4));
q_thr_4 = heat_flux(i_4);

pct_4 = 100*q_thr_4/max(heat_flux);

% Threshold 3
cnt_5 = 0;
thr_5 = 3*ten; % threshold
for i=1:size(moment_5')
    cnt_5 = cnt_5 + 1;
    if(moment_5(i)>= thr_5)
        break;
    end
end
t_thr_5 = tn(cnt_5);

[d_5, i_5] = min(abs(t - t_thr_5));
q_thr_5 = heat_flux(i_5);

pct_5 = 100*q_thr_5/max(heat_flux);

% Threshold 5
cnt_6 = 0;
thr_6 = 5*ten; % threshold
for i=1:size(moment_5')
    cnt_6 = cnt_6 + 1;
    if(moment_5(i)>= thr_6)
        break;
    end
end
t_thr_6 = tn(cnt_6);

[d_6, i_6] = min(abs(t - t_thr_6));
q_thr_6 = heat_flux(i_6);

pct_6 = 100*q_thr_6/max(heat_flux);



% Log Scale 
% Threshold 0.01
cnt_7 = 0;
thr_7 = 0.01*ten; % threshold
for i=1:size(moment_5')
    cnt_7 = cnt_7 + 1;
    if(moment_5(i)>= thr_7)
        break;
    end
end
t_thr_7 = tn(cnt_7);

[d_7, i_7] = min(abs(t - t_thr_7));
q_thr_7 = heat_flux(i_7);

pct_7 = 100*q_thr_7/max(heat_flux);

% Threshold 0.032
cnt_8 = 0;
thr_8 = 0.032*ten; % threshold
for i=1:size(moment_5')
    cnt_8 = cnt_8 + 1;
    if(moment_5(i)>= thr_8)
        break;
    end
end
t_thr_8 = tn(cnt_8);

[d_8, i_8] = min(abs(t - t_thr_8));
q_thr_8 = heat_flux(i_8);

pct_8 = 100*q_thr_8/max(heat_flux);


% Threshold 0.1
cnt_9 = 0;
thr_9 = 0.1*ten; % threshold
for i=1:size(moment_5')
    cnt_9 = cnt_9 + 1;
    if(moment_5(i)>= thr_9)
        break;
    end
end
t_thr_9 = tn(cnt_9);

[d_9, i_9] = min(abs(t - t_thr_9));
q_thr_9 = heat_flux(i_9);

pct_9 = 100*q_thr_9/max(heat_flux);

% Threshold 0.316
cnt_10 = 0;
thr_10 = 0.316*ten; % threshold
for i=1:size(moment_5')
    cnt_10 = cnt_10 + 1;
    if(moment_5(i)>= thr_10)
        break;
    end
end
t_thr_10 = tn(cnt_10);

[d_10, i_10] = min(abs(t - t_thr_10));
q_thr_10 = heat_flux(i_10);

pct_10 = 100*q_thr_10/max(heat_flux);


% Threshold 1
cnt_11 = 0;
thr_11 = 1*ten; % threshold
for i=1:size(moment_5')
    cnt_11 = cnt_11 + 1;
    if(moment_5(i)>= thr_11)
        break;
    end
end
t_thr_11 = tn(cnt_11);

[d_11, i_11] = min(abs(t - t_thr_11));
q_thr_11 = heat_flux(i_11);

pct_11 = 100*q_thr_11/max(heat_flux);

% Threshold 10
cnt_12 = 0;
thr_12 = 10*ten; % threshold
for i=1:size(moment_5')
    cnt_12 = cnt_12 + 1;
    if(moment_5(i)>= thr_12)
        break;
    end
end
t_thr_12 = tn(cnt_12);

[d_12, i_12] = min(abs(t - t_thr_12));
q_thr_12 = heat_flux(i_12);

pct_12 = 100*q_thr_12/max(heat_flux);

%% Writing to Excel 

% Create a cell array of variable names
varnames = {pct_1; pct_2; pct_3; pct_4; pct_5; pct_6; pct_7; pct_8; pct_9; pct_10; pct_11; pct_12};

% Write the variable names to the first row of the Excel file
xlswrite('data.xlsx', varnames, 'Sheet1', myrange);
%writematrix(varnames,'data.xlsx','Sheet', 1, 'Range', 'Q2')
%% PLotting 
figure(1);

plot(tn, moment_5)

hold on;
title('4th Moment vs time');
xlabel('t (s)');
ylabel('4th Moment ');
axis([min(tn) max(tn) min(moment_5) max(moment_5)])
plot(t_CHF*ones(size(tn)),moment_5', '-r') %CHF Line indicator
plot(t_thr_2*ones(size(tn)),moment_5', 'k') % line at which vibrations start 
hold off;


