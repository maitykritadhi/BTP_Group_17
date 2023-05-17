clear all;
clc;

filename = '250_4_formula.xls';
sheet = 1;

xlRange = 'A2:A35035';
ylRange = 'C2:C35035';
zlRange = 'D2:D35035';
x = xlsread(filename,sheet,xlRange);
y = xlsread(filename,sheet,ylRange);
z = xlsread(filename,sheet,zlRange);
for i = 2: length(x)
    dx = x(i)-x(i-1);
    dy = y(i)-y(i-1);
    dz = z(i)-z(i-1);
    d = dx*dx + dy*dy + dz*dz;
    d = sqrt(d);
    vec(i-1) = d;
end

t = linspace(0, 4380, 35033);

figure(1);

plot(t,vec);
xlabel('t (s)');
ylabel('delta distance F6 F59 F20');

hold on;

plot(2581*ones(size(t)),vec', '-g') %CHF Line indicator
plot(150*ones(size(t)),vec', '-r') % natural convection
% plot(2582*ones(size(t)),vec', '-g') %nucleate boiling
plot(3304*ones(size(t)),vec', '-r') %transition boiling
plot(4000*ones(size(t)),vec', '-r') %nucleate boiling
hold off;




% filename= 'F6_F53_F59.avi';  