%%
prompt1 = 'Enter a floating-point value of window size in milliseconds: ';
x = input(prompt1, 's');
x = str2double(x);
window_size = x;
x=x/1000;
prompt2 = 'Enter a floating-point value of window overlap in percentage: ';
y = input(prompt2, 's');
y = str2double(y);
window_overlap = y;
%%
fprintf('1-Maximum F1 \n')
fprintf('2-Mean F3 \n')
fprintf('3-Avg Signal Level F4 \n')
fprintf('4-Peak F5 \n')
fprintf('5-RMS F6 \n')
fprintf('6-Crest Factor F7 \n')
fprintf('7-Std Deviation F9 \n')
fprintf('8-Impulse factor F10 \n')
fprintf('9-CRIS F11 \n')
fprintf('10-Coeff of Variation F12 \n')
fprintf('11-Inverse Coeff of Variation F13 \n')
fprintf('12-Energy F14 \n')
fprintf('13-K-Factor F16 \n')
fprintf('14-Kurtosis F17 \n')
fprintf('15-Kurtosis Factor F18 \n')
fprintf('16-5th Moment F20 \n')
fprintf('17-6th Moment F21 \n')
fprintf('18-Square Mean Root F22 \n')
fprintf('19-Skewness F23 \n')
fprintf('20-Skewness Factor F24 \n')
fprintf('21-Shape Factor F25 \n')
fprintf('22-Mean Absolute Deviation F26 \n')
fprintf('23-Variance F27 \n')
fprintf('24-Margin Factor F29 \n')
fprintf('25-F30 \n')
fprintf('26-Histogram Upper Bound F31 \n')
fprintf('27-Histogram Lower Bound F32 \n')
fprintf('28-Clearance Factor F33 \n')
fprintf('29-Mobility F37 \n')
fprintf('30-Complexity F38 \n')
fprintf('31-LogLog ratio F51 \n')
fprintf('32-F52 \n')
fprintf('33-Minimum F53 \n')
fprintf('34-Most Frequent values F54 \n')
fprintf('35-F55 \n')
fprintf('36-Mean Absolute value F56 \n')
fprintf('37-F57 \n')
fprintf('38-F58 \n')
fprintf('39-Median value F59 \n')
%%
prompt3 = 'Enter a floating-point value of feature: ';
z = input(prompt3, 's');
z = str2double(z);
% [x, y, z] = input('Enter three values: ');
% x = str2double(x);
% y = str2double(y);
% z = str2double(z);

%%
[data,fs] = audioread("1-s2.0-S2666386421000722-mmc4.mp3");
% data- input data, fs-sampling data
% size(data)

N = length(data);
t = (0:N-1)/fs;

figure(1);

subplot(2,1,1);
plot(t,data');
title('AE Intensity vs time');
xlabel('t (s)');
ylabel('A (a.u.)');
ax = gca;
annotation('textbox',[ax.Position(1)-0.1,ax.Position(2)+ax.Position(4)-0.02,0.09,0.09],'String',['Window size: ' num2str(window_size) 'ms'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');
annotation('textbox',[ax.Position(1) + 0.6, ax.Position(2) + ax.Position(4) - 0.02, 0.09, 0.09], 'String', ['Window overlap: ' num2str(window_overlap) '%'], 'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left');

% subplot(3,1,2);
% plot(t,data')
% title('Plot 1');
% xlabel('x-axis label');
% ylabel('y-axis label');
% 
% annotation('textbox',[-0.1,0.1,0.1,0.1],'String',['Window size: ' num2str(window_size) 'ms'],'FitBoxToText','on');
% annotation('textbox',[0.25,0.25,0.1,0.1],'String',['Window overlap:' num2str(window_overlap) '%'],'FitBoxToText','on');

hold on;

plot(2581*ones(size(t)),data', '-g') %CHF Line indicator
plot(150*ones(size(t)),data', '-r') % natural convection
% plot(2582*ones(size(t)),data', '-g') %nucleate boiling
plot(3304*ones(size(t)),data', '-r') %transition boiling
plot(4000*ones(size(t)),data', '-r') %nucleate boiling
hold off;

%% 
% Set the window size (in seconds)
f_d = x;

% Convert the window size to samples
f_size = f_d * fs;

% Calculate the number of frames
n = length(data);
n_f = floor(n/f_size);
if (y==0)
    ovlp = 1;
else
    ovlp = round(100/(100-y));
end

% Extract the frames
temp = 0;
for i = 1 : ovlp*n_f
    frames(i,:) = data(temp + 1 : temp + f_size);
    temp = temp + floor(f_size/ovlp);
%     if temp>size(data)
%         break;
%     end
end  

[num_rows, num_columns] = size(frames);
tn = linspace(0,t(end),num_rows);
%%
% vector = frames(:,f_size/2);
% writematrix(vector,'X.xls','Sheet', 1, 'Range', 'A2');

%%
if z==1
    %F1
    max_frames = max(frames');
    % tn = (0:N-1)/fs;
    % MAX 1. F1
%     figure(2);
    subplot(2,1,2);
    plot(tn,max_frames)
    
    hold on;
    title('Maximum F1 vs time');
    xlabel('t (s)');
    ylabel('Maximum F1');
%     annotation('textbox',[0.25,0.2,0.1,0.1],'String','Feature name: Maximum F1','FitBoxToText','on');
    
    plot(2581*ones(size(tn)),max_frames', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),max_frames', '-r') % natural convection
    % plot(2582*ones(size(tn)),max_frames', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),max_frames', '-r') %transition boiling
    plot(4000*ones(size(tn)),max_frames', '-r') %nucleate boiling
    hold off;
    % writematrix(max_frames','X.xls','Sheet', 1, 'Range', 'B2')
    % xlabel('Time (s)')
    % ylabel('MAX(xi)')
    % title('Signal in Time Domain  max vs Time')
    % axis([0 4500 0 1.2])
elseif (z==2)
    %F3
    mean_frames = mean(frames');
%     figure(3);
    subplot(2,1,2);
    plot(tn,mean_frames)
    hold on;
    title('Mean F3 vs time');
    xlabel('t (s)');
    ylabel('Mean F3');
    axis([0 4500 min(mean_frames) max(mean_frames)])
    
    plot(2581*ones(size(tn)),mean_frames', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),mean_frames', '-r') % natural convection
    % plot(2582*ones(size(tn)),mean_frames', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),mean_frames', '-r') %transition boiling
    plot(4000*ones(size(tn)),mean_frames', '-r') %nucleate boiling
    hold off;
    writematrix(mean_frames','X.xls','Sheet', 1, 'Range', 'C2')
    % xlabel('Time (s)')
    % ylabel('MEAN(xi)')
    % title('Signal in Time Domain  mean vs Time')
elseif (z==3)
    %F4
    mean_frames = mean(frames');
    avgsignal_frames= mean_frames./(10^-6);
    avgsignal_frames = 20*log(avgsignal_frames);
    %ampl_frames = ampl_frames-40
%     figure(17);
    subplot(2,1,2);
    plot(tn, avgsignal_frames)
    hold on;
    title('Avg Signal Level F4 vs time');
    xlabel('t (s)');
    ylabel('Avg Signal Level F4');
    axis([0 4500 min(real(avgsignal_frames)) max(real(avgsignal_frames))])
    
    plot(2581*ones(size(tn)),avgsignal_frames', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),avgsignal_frames', '-r') % natural convection
    % plot(2582*ones(size(tn)),avgsignal_frames', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),avgsignal_frames', '-r') %transition boiling
    plot(4000*ones(size(tn)),avgsignal_frames', '-r') %nucleate boiling
    hold off;
    writematrix(avgsignal_frames','X.xls','Sheet', 1, 'Range', 'D2')
    % xlabel('Time (s)')
    % ylabel('F03')
    % title('Signal in Time Domain F04 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==4)
    % F5 Peak
    max_frames = max(frames');
    min_frames = min(frames');
    peak_cal_vec = (max_frames - min_frames)/2.;
%     figure(4);
    subplot(2,1,2);
    plot(tn,peak_cal_vec)
    hold on;
    title('Peak F5 vs time');
    xlabel('t (s)');
    ylabel('Peak F5 F4');
    axis([0 4500 min(peak_cal_vec) max(peak_cal_vec)])
    
    plot(2581*ones(size(tn)),peak_cal_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),peak_cal_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),peak_cal_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),peak_cal_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),peak_cal_vec', '-r') %nucleate boiling
    hold off;
    writematrix(peak_cal_vec','X.xls','Sheet', 1, 'Range', 'E2')
    % xlabel('Time (s)')
    % ylabel('PEAK(xi)')
    % title('Signal in Time Domain  Peak vs Time')
    axis([0 4500 0 1.2])
elseif (z==5)
    % RMS F6
    rms_vec = rms(frames');
%     figure(5);
    subplot(2,1,2);
    plot(tn,rms_vec)
    hold on;
    title('RMS F6 vs time');
    xlabel('t (s)');
    ylabel('RMS F6');
    axis([0 4500 min(rms_vec) max(rms_vec)])
    
    plot(2581*ones(size(tn)),rms_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),rms_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),rms_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),rms_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),rms_vec', '-r') %nucleate boiling
    hold off;
%     writematrix(rms_vec','X.xls','Sheet', 1, 'Range', 'F2')
    writematrix(rms_vec','250_4_formula.xls','Sheet', 1, 'Range', 'A2')
    % xlabel('Time (s)')
    % ylabel('RMS(xi)')
    % title('Signal in Time Domain  RMS vs Time')
elseif (z==6)
    % Crest Factor F7
    max_frames = max(frames');
    min_frames = min(frames');
    peak_cal_vec = (max_frames - min_frames)/2.;

    rms_vec = rms(frames');

    crest_fac_vec = peak_cal_vec./rms_vec;
%     figure(6);
    subplot(2,1,2);
    plot(tn,crest_fac_vec)
    hold on;
    title('Crest Factor F7 vs time');
    xlabel('t (s)');
    ylabel('Crest Factor F7');
    axis([0 4500 min(crest_fac_vec) max(crest_fac_vec)])
    
    plot(2581*ones(size(tn)),crest_fac_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),crest_fac_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),crest_fac_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),crest_fac_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),crest_fac_vec', '-r') %nucleate boiling
    hold off;
    writematrix(crest_fac_vec','X.xls','Sheet', 1, 'Range', 'G2')
    % xlabel('Time (s)')
    % ylabel('CREST_FACTOR(xi)')
    % title('Signal in Time Domain  CREST FACTOR vs Time')
elseif (z==7)
    %STANDARD DEVIATION F9
    std_vec = std(frames');
%     figure(7);
    subplot(2,1,2);
    plot(tn,std_vec)
    hold on;
    title('Std Deviation F9 vs time');
    xlabel('t (s)');
    ylabel('Std Deviation F9');
    axis([0 4500 min(std_vec) max(std_vec)])
    
    plot(2581*ones(size(tn)),std_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),std_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),std_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),std_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),std_vec', '-r') %nucleate boiling
    hold off;
    writematrix(std_vec','X.xls','Sheet', 1, 'Range', 'H2')
    % xlabel('Time (s)')
    % ylabel('STANDARD_DEVIATION(xi)')
    % title('Signal in Time Domain  STANDARD DEVIATION vs Time')
elseif (z==8)
    %F10(Impulse Factor)
    max_frames = max(frames');
    mean_frames = mean(frames');

    impulse = max_frames./mean_frames;
%     figure(18);
    subplot(2,1,2);
    plot(tn, impulse)
    hold on;
    title('Impulse factor F10 vs time');
    xlabel('t (s)');
    ylabel('Impulse factor F10');
    axis([0 4500 min(impulse) max(impulse)])
    
    plot(2581*ones(size(tn)),impulse', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),impulse', '-r') % natural convection
    % plot(2582*ones(size(tn)),impulse', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),impulse', '-r') %transition boiling
    plot(4000*ones(size(tn)),impulse', '-r') %nucleate boiling
    hold off;
    writematrix(impulse','X.xls','Sheet', 1, 'Range', 'I2')
    % xlabel('Time (s)')
    % ylabel('F10')
    % title('Signal in Time Domain Impulse factor vs Time')
    % axis([0 1197 0 0.25])
elseif (z==9)
    %F11 CRIS
    max_frames = max(frames');
    min_frames = min(frames');
    peak_cal_vec = (max_frames - min_frames)/2.;
    rms_vec = rms(frames');
    crest_fac_vec = peak_cal_vec./rms_vec;

    mean_frames = mean(frames');
    impulse = max_frames./mean_frames;

    std_vec = std(frames');

    cris_a = crest_fac_vec./(crest_fac_vec-2);
    cris_b = rms_vec.^2;
    id = ones(1,n_f*ovlp);
    cris_b = id./cris_b;
    cris_b = cris_b + impulse;
    cris_b = cris_b .* std_vec;
    cris = cris_a.*cris_b;
    cris = sqrt(cris);
    cris = (10.*rms_vec) + cris;
%     figure(19);
    subplot(2,1,2);
    plot(tn, cris)
    hold on;
    title('CRIS F11 vs time');
    xlabel('t (s)');
    ylabel('CRIS F11');
    axis([0 4500 min(cris) max(cris)])
    
    plot(2581*ones(size(tn)),cris', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),cris', '-r') % natural convection
    % plot(2582*ones(size(tn)),cris', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),cris', '-r') %transition boiling
    plot(4000*ones(size(tn)),cris', '-r') %nucleate boiling
    hold off;
    writematrix(cris','X.xls','Sheet', 1, 'Range', 'J2')
    % xlabel('Time (s)')
    % ylabel('F11')
    % title('Signal in Time Domain CRIS vs Time')
    % axis([0 1197 0 0.25])
elseif (z==10)
    %F12 coeff of var
    std_vec = std(frames');
    mean_frames = mean(frames');

    coeff_var = std_vec.^2;
    coeff_var = coeff_var./mean_frames;
%     figure(20);
    subplot(2,1,2);
    plot(tn, coeff_var)
    hold on;
    title('Coeff of Variation F12 vs time');
    xlabel('t (s)');
    ylabel('Coeff of Variation F12');
    axis([0 4500 min(coeff_var) max(coeff_var)])
    
    plot(2581*ones(size(tn)),coeff_var', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),coeff_var', '-r') % natural convection
    % plot(2582*ones(size(tn)),coeff_var', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),coeff_var', '-r') %transition boiling
    plot(4000*ones(size(tn)),coeff_var', '-r') %nucleate boiling
    hold off;
    writematrix(coeff_var','X.xls','Sheet', 1, 'Range', 'K2')
    % xlabel('Time (s)')
    % ylabel('F12')
    % title('Signal in Time Domain F12 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==11)
    %F13 Inverse coeff of var
    std_vec = std(frames');
    mean_frames = mean(frames');
    coeff_var = std_vec.^2;
    coeff_var = coeff_var./mean_frames;

    inv_coeff_var = 1./coeff_var;
%     figure(21);
    subplot(2,1,2);
    plot(tn, inv_coeff_var)
    hold on;
    title('Inverse Coeff of Variation F13 vs time');
    xlabel('t (s)');
    ylabel('Inverse Coeff of Variation F13');
    axis([0 4500 min(inv_coeff_var) max(inv_coeff_var)])
    
    plot(2581*ones(size(tn)),inv_coeff_var', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),inv_coeff_var', '-r') % natural convection
    % plot(2582*ones(size(tn)),inv_coeff_var', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),inv_coeff_var', '-r') %transition boiling
    plot(4000*ones(size(tn)),inv_coeff_var', '-r') %nucleate boiling
    hold off;
    writematrix(inv_coeff_var','X.xls','Sheet', 1, 'Range', 'L2')
    % xlabel('Time (s)')
    % ylabel('F`13')
    % title('Signal in Time Domain F13 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==12)
    %Energy F14
    f_2 = frames.^2;
    energy_vec = sum(f_2');
%     figure(10);
    subplot(2,1,2);
    plot(tn,energy_vec)
    hold on;
    title('Energy F14 vs time');
    xlabel('t (s)');
    ylabel('Energy F14');
    axis([0 4500 min(energy_vec) max(energy_vec)])
    
    plot(2581*ones(size(tn)),energy_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),energy_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),energy_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),energy_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),energy_vec', '-r') %nucleate boiling
    hold off;
    writematrix(energy_vec','X.xls','Sheet', 1, 'Range', 'M2')
    % xlabel('Time (s)')
    % ylabel('Energy')
    % title('Signal in Time Domain  Energy vs Time')
elseif (z==13)
    % K-Factor F16
    max_frames = max(frames');
    min_frames = min(frames');
    peak_cal_vec = (max_frames - min_frames)/2.;
    rms_vec = rms(frames');

    kfactor_vec = peak_cal_vec.*rms_vec;
%     figure(8);
    subplot(2,1,2);
    plot(tn,kfactor_vec)
    hold on;
    title('K-Factor F16 vs time');
    xlabel('t (s)');
    ylabel('K-Factor F16');
    axis([0 4500 min(kfactor_vec) max(kfactor_vec)])
    
    plot(2581*ones(size(tn)),std_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),std_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),std_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),std_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),std_vec', '-r') %nucleate boiling
    hold off;
    writematrix(std_vec','X.xls','Sheet', 1, 'Range', 'N2')
    % xlabel('Time (s)')
    % ylabel('K-Factor')
    % title('Signal in Time Domain  KFactor vs Time')
elseif (z==14)
    %Kurtosis F17
    mean_frames = mean(frames');
    std_vec = std(frames');

    sub_kurt = frames - mean_frames';
    pow_kurt = sub_kurt.^4;
    sum_kurt = sum(pow_kurt');
    sum_kurt1 = sum(pow_kurt')./f_size;
    f9_4 = std_vec.^4;
    kurtosis = sum_kurt1./f9_4;
%     figure(9);
    subplot(2,1,2);
    plot(tn,kurtosis)
    hold on;
    title('Kurtosis F17 vs time');
    xlabel('t (s)');
    ylabel('Kurtosis F17');
    axis([0 4500 min(kurtosis) max(kurtosis)])
    
    plot(2581*ones(size(tn)),kurtosis', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),kurtosis', '-r') % natural convection
    % plot(2582*ones(size(tn)),kurtosis', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),kurtosis', '-r') %transition boiling
    plot(4000*ones(size(tn)),kurtosis', '-r') %nucleate boiling
    hold off;
    writematrix(kurtosis','X.xls','Sheet', 1, 'Range', 'O2')
    % xlabel('Time (s)')
    % ylabel('Kurtosis')
    % title('Signal in Time Domain Kurtosis vs Time')
elseif (z==15)
    % F18 Kurtosis Factor
    rms_vec = rms(frames');

    mean_frames = mean(frames');
    std_vec = std(frames');
    sub_kurt = frames - mean_frames';
    pow_kurt = sub_kurt.^4;
    sum_kurt = sum(pow_kurt');
    sum_kurt1 = sum(pow_kurt')./f_size;
    f9_4 = std_vec.^4;
    kurtosis = sum_kurt1./f9_4;

    kurt_fact = kurtosis./(rms_vec.^4);
%     figure(22);
    subplot(2,1,2);
    plot(tn, kurt_fact)
    hold on;
    title('Kurtosis Factor F18 vs time');
    xlabel('t (s)');
    ylabel('Kurtosis Factor F18');
    axis([0 4500 min(kurt_fact) max(kurt_fact)])
    
    plot(2581*ones(size(tn)),kurt_fact', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),kurt_fact', '-r') % natural convection
    % plot(2582*ones(size(tn)),kurt_fact', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),kurt_fact', '-r') %transition boiling
    plot(4000*ones(size(tn)),kurt_fact', '-r') %nucleate boiling
    hold off;
    writematrix(kurt_fact','X.xls','Sheet', 1, 'Range', 'P2')
    writematrix(kurt_fact','250_4_formula.xls','Sheet', 1, 'Range', 'F2')
    % xlabel('Time (s)')
    % ylabel('F18')
    % title('Signal in Time Domain F18 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==16)
    % F20 
    %frames1= frames - 
    mean_frames = mean(frames');


    for i = 1:n_f*ovlp
        frames1(i,:) = frames(i,:)-mean_frames(i);
    end
    frames1 = frames1.^4;
    %col = size(frames1);
    %col = 110250;
    [rows,col] = size(frames1);
    
    moment_5 = sum(frames1');
%     figure(23);
    subplot(2,1,2);
    plot(tn, moment_5)
    
    hold on;
    title('5th Moment F20 vs time');
    xlabel('t (s)');
    ylabel('5th Moment F20');
    axis([0 4500 min(moment_5) max(moment_5)])
    
    plot(2581*ones(size(tn)),moment_5', '-r','Linewidth',1.1) %CHF Line indicator
%     plot(150*ones(size(tn)),moment_5', '-r') % natural convection
%     % plot(2582*ones(size(tn)),moment_5', '-g') %nucleate boiling
%     plot(3304*ones(size(tn)),moment_5', '-r') %transition boiling
%     plot(4000*ones(size(tn)),moment_5', '-r') %nucleate boiling
    hold off;
    %writematrix(moment_5','X.xls','Sheet', 1, 'Range', 'Q2')
    writematrix(moment_5','250_4_formula.xls','Sheet', 1, 'Range', 'D2')
    % xlabel('Time (s)')
    % ylabel('F20')
    % title('Signal in Time Domain F20 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==17)
    % F21
    %frames1= frames - 
    mean_frames = mean(frames');

    for i = 1:n_f*ovlp
        frames2(i,:) = frames(i,:)-mean_frames(i);
    end
    frames2 = frames2.^6;
    %col = size(frames1);
    %col = 110250;
    %[rows,col] = size(frames1);
    
    moment_6 = sum(frames2');
%     figure(24);
    subplot(2,1,2);
    plot(tn, moment_6)
    hold on;
    title('6th Moment F21 vs time');
    xlabel('t (s)');
    ylabel('6th Moment F21');
    axis([0 4500 min(moment_6) max(moment_6)])
    
    plot(2581*ones(size(tn)),moment_6', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),moment_6', '-r') % natural convection
    % plot(2582*ones(size(tn)),moment_6', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),moment_6', '-r') %transition boiling
    plot(4000*ones(size(tn)),moment_6', '-r') %nucleate boiling
    hold off;
    writematrix(moment_6','X.xls','Sheet', 1, 'Range', 'R2')
    % xlabel('Time (s)')
    % ylabel('F21')
    % title('Signal in Time Domain F21 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==18)
    %Square Mean root F22
    abs_frames = abs(frames);
    sq_vec = abs_frames.^0.5;  % taken absolute values
    add_vec_by_N = sum(sq_vec')./f_size;
    smr_vec = add_vec_by_N.^2;
%     figure(11);
    subplot(2,1,2);
    plot(tn,smr_vec)
    hold on;
    title('Square Mean Root F22 vs time');
    xlabel('t (s)');
    ylabel('Square Mean Root F22');
    axis([0 4500 min(smr_vec) max(smr_vec)])
    
    plot(2581*ones(size(tn)),smr_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),smr_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),smr_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),smr_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),smr_vec', '-r') %nucleate boiling
    hold off;
    writematrix(smr_vec','X.xls','Sheet', 1, 'Range', 'S2')
    % xlabel('Time (s)')
    % ylabel('Square Mean Root')
    % title('Signal in Time Domain Square Mean Root vs Time')
    % axis([0 1197 0 0.16])
elseif (z==19)
    % Skewness F23
    mean_frames = mean(frames');
    std_vec = std(frames');
    sub_kurt = frames - mean_frames';

    pow_kurt1 = sub_kurt.^3;
    sum_kurt0 = sum(pow_kurt1');
    sum_kurt0 = sum_kurt0./f_size;
    f9_3 = std_vec.^3;
    skewness = sum_kurt0./f9_3;
%     figure(25);
    subplot(2,1,2);
    plot(tn,skewness)
    hold on;
    title('Skewness F23 vs time');
    xlabel('t (s)');
    ylabel('Skewness F23');
    axis([0 4500 min(skewness) max(skewness)])
    
    plot(2581*ones(size(tn)),skewness', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),skewness', '-r') % natural convection
    % plot(2582*ones(size(tn)),skewness', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),skewness', '-r') %transition boiling
    plot(4000*ones(size(tn)),skewness', '-r') %nucleate boiling
    hold off;
    writematrix(skewness','X.xls','Sheet', 1, 'Range', 'T2')
    % xlabel('Time (s)')
    % ylabel('Skewness')
    % title('Signal in Time Domain Skewness vs Time')
    % axis([0 1197 0 0.25])
elseif (z==20)
    % Skewness Factor F24
    rms_vec = rms(frames');

    mean_frames = mean(frames');
    std_vec = std(frames');
    sub_kurt = frames - mean_frames';
    pow_kurt1 = sub_kurt.^3;
    sum_kurt0 = sum(pow_kurt1');
    sum_kurt0 = sum_kurt0./f_size;
    f9_3 = std_vec.^3;
    skewness = sum_kurt0./f9_3;

    skewness_factor = skewness./(rms_vec.^3);
%     figure(26);
    subplot(2,1,2);
    plot(tn,skewness_factor)
    hold on;
    title('Skewness Factor F24 vs time');
    xlabel('t (s)');
    ylabel('Skewness Factor F24');
    axis([0 4500 min(skewness_factor) max(skewness_factor)])
    
    plot(2581*ones(size(tn)),skewness_factor', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),skewness_factor', '-r') % natural convection
    % plot(2582*ones(size(tn)),skewness_factor', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),skewness_factor', '-r') %transition boiling
    plot(4000*ones(size(tn)),skewness_factor', '-r') %nucleate boiling
    hold off;
    writematrix(skewness_factor','X.xls','Sheet', 1, 'Range', 'U2')
    % xlabel('Time (s)')
    % ylabel('Skewness Factor')
    % title('Signal in Time Domain Skewness Factor vs Time')
    % axis([0 1197 0 0.25])
elseif (z==21)
    % F25 Shape Factor
    rms_vec = rms(frames');
    mean_frames = mean(frames');

    shape_factor = rms_vec./mean_frames;
%     figure(27);
    subplot(2,1,2);
    plot(tn,shape_factor)
    hold on;
    title('Shape Factor F25 vs time');
    xlabel('t (s)');
    ylabel('Shape Factor F25');
    axis([0 4500 min(shape_factor) max(shape_factor)])
    
    plot(2581*ones(size(tn)),shape_factor', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),shape_factor', '-r') % natural convection
    % plot(2582*ones(size(tn)),shape_factor', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),shape_factor', '-r') %transition boiling
    plot(4000*ones(size(tn)),shape_factor', '-r') %nucleate boiling
    hold off;
    writematrix(shape_factor','X.xls','Sheet', 1, 'Range', 'V2')
    % xlabel('Time (s)')
    % ylabel('Shape Factor')
    % title('Signal in Time Domain Shape Factor vs Time')
    % axis([0 1197 0 0.25])
elseif (z==22)
    % F26 Mean Abs Dev
    mean_frames = mean(frames');

    mean_abs = frames - mean_frames';
    mean_abs = abs(mean_abs);
    mean_abs = sum(mean_abs');
    mean_abs = mean_abs./f_size;
%     figure(28);
    subplot(2,1,2);
    plot(tn,mean_abs)
    hold on;
    title('Mean Absolute Deviation F26 vs time');
    xlabel('t (s)');
    ylabel('Mean Absolute Deviation F26');
    axis([0 4500 min(mean_abs) max(mean_abs)])
    
    plot(2581*ones(size(tn)),mean_abs', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),mean_abs', '-r') % natural convection
    % plot(2582*ones(size(tn)),mean_abs', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),mean_abs', '-r') %transition boiling
    plot(4000*ones(size(tn)),mean_abs', '-r') %nucleate boiling
    hold off;
    writematrix(mean_abs','X.xls','Sheet', 1, 'Range', 'W2')
    % xlabel('Time (s)')
    % ylabel('Mean Abs Dev')
    % title('Signal in Time Domain Mean Abs Dev vs Time')
    % axis([0 1197 0 0.25])
elseif (z==23)
    % F27 Variance
    mean_frames = mean(frames');
    sub_kurt = frames - mean_frames';

    variance = sub_kurt.^2;
    variance = sum(variance');
    variance = variance./f_size;
%     figure(29);
    subplot(2,1,2);
    plot(tn,variance)
    hold on;
    title('Variance F27 vs time');
    xlabel('t (s)');
    ylabel('Variance F27');
    axis([0 4500 min(variance) max(variance)])
    
    plot(2581*ones(size(tn)),variance', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),variance', '-r') % natural convection
    % plot(2582*ones(size(tn)),variance', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),variance', '-r') %transition boiling
    plot(4000*ones(size(tn)),variance', '-r') %nucleate boiling
    hold off;
    writematrix(variance','X.xls','Sheet', 1, 'Range', 'X2')
    % xlabel('Time (s)')
    % ylabel('Variance,F27')
    % title('Signal in Time Domain Variance vs Time')
    % axis([0 1197 0 0.25])
elseif (z==24)
    %F28/F29 
    %FM4 = F1/F22  
    abs_frames = abs(frames);
    sq_vec = abs_frames.^0.5;  % taken absolute values
    add_vec_by_N = sum(sq_vec')./f_size;
    smr_vec = add_vec_by_N.^2;

    max_frames = max(frames');

    fm4_vec = max_frames./smr_vec;
%     figure(12);
    subplot(2,1,2);
    plot(tn,fm4_vec)
    hold on;
    title('Margin Factor F29 vs time');
    xlabel('t (s)');
    ylabel('Margin Factor F29');
    axis([0 4500 min(fm4_vec) max(fm4_vec)])
    
    plot(2581*ones(size(tn)),fm4_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),fm4_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),fm4_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),fm4_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),fm4_vec', '-r') %nucleate boiling
    hold off;
    writematrix(fm4_vec','X.xls','Sheet', 1, 'Range', 'Y2')
    % xlabel('Time (s)')
    % ylabel('F28/F29')
    % title('Signal in Time Domain F28/F29 vs Time')
    % axis([0 1197 0 4000])
elseif (z==25)
    % F30 
    abs_frames = abs(frames);
    sq_vec = abs_frames.^0.5;  % taken absolute values
    add_vec_by_N = sum(sq_vec')./f_size;
    smr_vec = add_vec_by_N.^2;

    std_vec = std(frames');

    f_30 = smr_vec./std_vec;
    f_30 = sqrt(f_30);
%     figure(30);
    subplot(2,1,2);
    plot(tn,f_30)
    hold on;
    title('F30 vs time');
    xlabel('t (s)');
    ylabel('F30');
    axis([0 4500 min(f_30) max(f_30)])
    
    plot(2581*ones(size(tn)),f_30', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),f_30', '-r') % natural convection
    % plot(2582*ones(size(tn)),f_30', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),f_30', '-r') %transition boiling
    plot(4000*ones(size(tn)),f_30', '-r') %nucleate boiling
    hold off;
    writematrix(f_30','X.xls','Sheet', 1, 'Range', 'Z2')
    % xlabel('Time (s)')
    % ylabel('F 30')
    % title('Signal in Time Domain F 30vs Time')
    % axis([0 1197 0 0.25])
elseif (z==26)
    % F31
    % Upper Bound 
    max_frames = max(frames');
    min_frames = min(frames');

    upper_bound = max_frames-min_frames;
    upper_bound = upper_bound./(2*(f_size-1));
    upper_bound = upper_bound +max_frames;
%     figure(31);
    subplot(2,1,2);
    plot(tn,upper_bound)
    hold on;
    title('Histogram Upper Bound F31 vs time');
    xlabel('t (s)');
    ylabel('Histogram Upper Bound F31');
    axis([0 4500 min(upper_bound) max(upper_bound)])
    
    plot(2581*ones(size(tn)),upper_bound', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),upper_bound', '-r') % natural convection
    % plot(2582*ones(size(tn)),upper_bound', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),upper_bound', '-r') %transition boiling
    plot(4000*ones(size(tn)),upper_bound', '-r') %nucleate boiling
    hold off;
    writematrix(upper_bound','X.xls','Sheet', 1, 'Range', 'AA2')
    % xlabel('Time (s)')
    % ylabel('Histogram Upper Bound')
    % title('Signal in Time Domain Histogram Upper Bound vs Time')
    % axis([0 1197 0 0.25])
elseif (z==27)
    %F32
    % lower bound
    max_frames = max(frames');
    min_frames = min(frames');

    lower_bound = max_frames-min_frames;
    lower_bound = lower_bound./(2*(f_size-1));
    lower_bound = min_frames - lower_bound;
%     figure(32);
    subplot(2,1,2);
    plot(tn,lower_bound)
    hold on;
    title('Histogram Lower Bound F32 vs time');
    xlabel('t (s)');
    ylabel('Histogram Lower Bound F32');
    axis([0 4500 min(lower_bound) max(lower_bound)])
    
    plot(2581*ones(size(tn)),lower_bound', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),lower_bound', '-r') % natural convection
    % plot(2582*ones(size(tn)),lower_bound', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),lower_bound', '-r') %transition boiling
    plot(4000*ones(size(tn)),lower_bound', '-r') %nucleate boiling
    hold off;
    writematrix(lower_bound','X.xls','Sheet', 1, 'Range', 'AB2')
    % xlabel('Time (s)')
    % ylabel('Histogram Lower Bound')
    % title('Signal in Time Domain Histogram Lower Bound vs Time')
    % axis([0 1197 0 0.25])
elseif (z==28)
    % F33 Clearnace Factor 
    max_frames = max(frames');
    min_frames = min(frames');
    peak_cal_vec = (max_frames - min_frames)/2.;

    clearance = abs(frames);
    clearance = sum(clearance');
    clearance = clearance./f_size;
    clearance = clearance.^2;
    clearance = peak_cal_vec./clearance;
%     figure(33);
    subplot(2,1,2);
    plot(tn,clearance)
    hold on;
    title('Clearance Factor F33 vs time');
    xlabel('t (s)');
    ylabel('Clearance Factor F33');
    axis([0 4500 min(clearance) max(clearance)])
    
    plot(2581*ones(size(tn)),clearance', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),clearance', '-r') % natural convection
    % plot(2582*ones(size(tn)),clearance', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),clearance', '-r') %transition boiling
    plot(4000*ones(size(tn)),clearance', '-r') %nucleate boiling
    hold off;
    writematrix(clearance','X.xls','Sheet', 1, 'Range', 'AC2')
    % xlabel('Time (s)')
    % ylabel('Clearance')
    % title('Signal in Time Domain Clearance vs Time')
    % axis([0 1197 0 0.25])
elseif (z==29)
    % F37
    std_vec = std(frames');

    derivative_frames = diff(frames);
    std_derivative_frames = std(derivative_frames');
    std_derivative_frames(end+1) = 0;
    ans = std_derivative_frames/std_vec;
    std_derivative_frames = std_derivative_frames.*ans;
%     figure(34);
    subplot(2,1,2);
    plot(tn,std_derivative_frames)
    hold on;
    title('Mobility F37 vs time');
    xlabel('t (s)');
    ylabel('Mobility F37');
    axis([0 4500 min(std_derivative_frames) max(std_derivative_frames)])
    
    plot(2581*ones(size(tn)),std_derivative_frames', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),std_derivative_frames', '-r') % natural convection
    % plot(2582*ones(size(tn)),std_derivative_frames', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),std_derivative_frames', '-r') %transition boiling
    plot(4000*ones(size(tn)),std_derivative_frames', '-r') %nucleate boiling
    hold off;
    writematrix(std_derivative_frames','X.xls','Sheet', 1, 'Range', 'AD2')
    % xlabel('Time (s)')
    % ylabel('F 34')
    % title('Signal in F 34 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==30)
    % F38
    std_vec = std(frames');
    derivative_frames = diff(frames);
    std_derivative_frames = std(derivative_frames');
    std_derivative_frames(end+1) = 0;
    ans = std_derivative_frames/std_vec;
    std_derivative_frames = std_derivative_frames.*ans;

    double_derivative = diff(derivative_frames);
    std_double_derivative = std(double_derivative');
    std_double_derivative(end+1) = 0;
    std_double_derivative(end+1) = 0;
    ans_1 = std_double_derivative/std_derivative_frames;
    ans_2 = std_derivative_frames/std_vec;
    f38 = ans_1/ans_2;
    f38 = f38*std_double_derivative;
%     figure(35);
    subplot(2,1,2);
    plot(tn,f38)
    hold on;
    title('Complexity F38 vs time');
    xlabel('t (s)');
    ylabel('Complexity F38');
    axis([0 4500 min(f38) max(f38)])
    
    plot(2581*ones(size(tn)),f38', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),f38', '-r') % natural convection
    % plot(2582*ones(size(tn)),f38', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),f38', '-r') %transition boiling
    plot(4000*ones(size(tn)),f38', '-r') %nucleate boiling
    hold off;
    writematrix(f38','X.xls','Sheet', 1, 'Range', 'AE2')
    % xlabel('Time (s)')
    % ylabel('F 38')
    % title('Signal in F 38 vs Time')
    % axis([0 1197 0 0.25])
elseif(z==31)
    % F51 Log-Log Ratio
    std_vec = std(frames');

    log_frames = abs(frames);
    log_frames = frames + 1;
    log_frames = log(log_frames);
    log_frames = sum(log_frames');
    log_a = log(std_vec);
    log_frames = log_frames./log_a;
%     figure(43);
    subplot(2,1,2);
    plot(tn, log_frames)
    hold on;
    title('LogLog ratio F51 vs time');
    xlabel('t (s)');
    ylabel('LogLog ratio F51');
    axis([0 4500 min(log_frames) max(log_frames)])
    
    plot(2581*ones(size(tn)),log_frames', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),log_frames', '-r') % natural convection
    % plot(2582*ones(size(tn)),log_frames', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),log_frames', '-r') %transition boiling
    plot(4000*ones(size(tn)),log_frames', '-r') %nucleate boiling
    hold off;
    writematrix(log_frames','X.xls','Sheet', 1, 'Range', 'AF2')
    % xlabel('Time (s)')
    % ylabel('F 41')
    % title('Signal in F 41 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==32)
    % F52
    rms_vec = rms(frames');

    mean_frames = mean(frames');
    std_vec = std(frames');
    sub_kurt = frames - mean_frames';
    pow_kurt = sub_kurt.^4;
    sum_kurt = sum(pow_kurt');
    sum_kurt1 = sum(pow_kurt')./f_size;
    f9_4 = std_vec.^4;
    kurtosis = sum_kurt1./f9_4;

    f52 = rms_vec.*kurtosis;
%     figure(42);
    subplot(2,1,2);
    plot(tn,f52)
    hold on;
    title('F52 vs time');
    xlabel('t (s)');
    ylabel('F52');
    axis([0 4500 min(f52) max(f52)])
    
    plot(2581*ones(size(tn)),f52', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),f52', '-r') % natural convection
    % plot(2582*ones(size(tn)),f52', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),f52', '-r') %transition boiling
    plot(4000*ones(size(tn)),f52', '-r') %nucleate boiling
    hold off;
    writematrix(f52','X.xls','Sheet', 1, 'Range', 'AG2')
    % xlabel('Time (s)')
    % ylabel('F 41')
    % title('Signal in F 41 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==33)
    % F53
    min_frames = min(frames');

%     figure(36);
    subplot(2,1,2);
    plot(tn, min_frames)
    hold on;
    title('F53 vs time');
    xlabel('t (s)');
    ylabel('F53');
    axis([0 4500 min(min_frames) max(min_frames)])
    
    plot(2581*ones(size(tn)),min_frames', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),min_frames', '-r') % natural convection
    % plot(2582*ones(size(tn)),min_frames', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),min_frames', '-r') %transition boiling
    plot(4000*ones(size(tn)),min_frames', '-r') %nucleate boiling
    hold off;
    %writematrix(min_frames','X.xls','Sheet', 1, 'Range', 'AH2')
    writematrix(min_frames','250_4_formula.xls','Sheet', 1, 'Range', 'B2')
    % xlabel('Time (s)')
    % ylabel('F 53')
    % title('Signal in F 53 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==34)
    % F54
    mode_frames = mode(frames');
%     figure(37);
    subplot(2,1,2);
    plot(tn,mode_frames)
    hold on;
    title('F54 vs time');
    xlabel('t (s)');
    ylabel('F54');
    axis([0 4500 min(mode_frames) max(mode_frames)])
    
    plot(2581*ones(size(tn)),mode_frames', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),mode_frames', '-r') % natural convection
    % plot(2582*ones(size(tn)),mode_frames', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),mode_frames', '-r') %transition boiling
    plot(4000*ones(size(tn)),mode_frames', '-r') %nucleate boiling
    hold off;
    writematrix(mode_frames','X.xls','Sheet', 1, 'Range', 'AI2')
    % xlabel('Time (s)')
    % ylabel('F 41')
    % title('Signal in F 41 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==35)
    % F55
    rms_vec = rms(frames');
    max_frames = max(frames');

    f55 = max_frames./rms_vec;
%     figure(38);
    subplot(2,1,2);
    plot(tn,f55)
    hold on;
    title('F55 vs time');
    xlabel('t (s)');
    ylabel('F55');
    axis([0 4500 min(f55) max(f55)])
    
    plot(2581*ones(size(tn)),f55', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),f55', '-r') % natural convection
    % plot(2582*ones(size(tn)),f55', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),f55', '-r') %transition boiling
    plot(4000*ones(size(tn)),f55', '-r') %nucleate boiling
    hold off;
    writematrix(f55','X.xls','Sheet', 1, 'Range', 'AJ2')
    % xlabel('Time (s)')
    % ylabel('F 41')
    % title('Signal in F 41 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==36)
    % Mean Absolute Value F56
    abs_vec = abs(frames);
    sum_abs_vec = sum(abs_vec')./f_size;
%     figure(13);
    subplot(2,1,2);
    plot(tn,sum_abs_vec)
    hold on;
    title('Mean Absolute value F56 vs time');
    xlabel('t (s)');
    ylabel('Mean Absolute value F56');
    axis([0 4500 min(sum_abs_vec) max(sum_abs_vec)])
    
    plot(2581*ones(size(tn)),sum_abs_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),sum_abs_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),sum_abs_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),sum_abs_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),sum_abs_vec', '-r') %nucleate boiling
    hold off;
    writematrix(sum_abs_vec','X.xls','Sheet', 1, 'Range', 'AK2')
    % xlabel('Time (s)')
    % ylabel('Mean Absolute Value')
    % title('Signal in Time Domain Mean Absolute Value vs Time')
    % axis([0 1197 0 0.25])
elseif (z==37)
    % F57
    rms_vec = rms(frames');
    abs_vec = abs(frames);
    sum_abs_vec = sum(abs_vec')./f_size;

    f57 = sum_abs_vec./rms_vec;
%     figure(39);
    subplot(2,1,2);
    plot(tn, f57)
    hold on;
    title('F57 vs time');
    xlabel('t (s)');
    ylabel('F57');
    axis([0 4500 min(f57) max(f57)])
    
    plot(2581*ones(size(tn)),f57', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),f57', '-r') % natural convection
    % plot(2582*ones(size(tn)),f57', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),f57', '-r') %transition boiling
    plot(4000*ones(size(tn)),f57', '-r') %nucleate boiling
    hold off;
    writematrix(f57','X.xls','Sheet', 1, 'Range', 'AL2')
    % xlabel('Time (s)')
    % ylabel('F 41')
    % title('Signal in F 41 vs Time')
    % axis([0 1197 0 0.25])
elseif (z==38)
    % F58
    max_frames = max(frames');
    abs_vec = abs(frames);
    sum_abs_vec = sum(abs_vec')./f_size;

    f_58_vec = sum_abs_vec./max_frames;
%     figure(15);
    subplot(2,1,2);
    plot(tn,f_58_vec)
    hold on;
    title('F58 vs time');
    xlabel('t (s)');
    ylabel('F58');
    axis([0 4500 min(f_58_vec) max(f_58_vec)])
    
    plot(2581*ones(size(tn)),f_58_vec', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),f_58_vec', '-r') % natural convection
    % plot(2582*ones(size(tn)),f_58_vec', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),f_58_vec', '-r') %transition boiling
    plot(4000*ones(size(tn)),f_58_vec', '-r') %nucleate boiling
    hold off;
    writematrix(f_58_vec','X.xls','Sheet', 1, 'Range', 'AM2')
    % xlabel('Time (s)')
    % ylabel('F58')
    % title('Signal in Time Domain F58 vs Time')
    % axis([0 1197 0 0.25])
else
    % F59 
    median_frames = median(frames');
%     figure(41);
    subplot(2,1,2);
    plot(tn,median_frames)
    hold on;
    title('Median value F59 vs time');
    xlabel('t (s)');
    ylabel('Median value F59');
    axis([0 4500 min(median_frames) max(median_frames)])
    
    plot(2581*ones(size(tn)),median_frames', '-g') %CHF Line indicator
    plot(150*ones(size(tn)),median_frames', '-r') % natural convection
    % plot(2582*ones(size(tn)),median_frames', '-g') %nucleate boiling
    plot(3304*ones(size(tn)),median_frames', '-r') %transition boiling
    plot(4000*ones(size(tn)),median_frames', '-r') %nucleate boiling
    hold off;
    writematrix(median_frames','X.xls','Sheet', 1, 'Range', 'AN2')
    writematrix(median_frames','250_4_formula.xls','Sheet', 1, 'Range', 'C2')
    % xlabel('Time (s)')
    % ylabel('F 59')
    % title('Signal in F 59 vs Time')
    % axis([0 1197 0 0.25])
end

    





