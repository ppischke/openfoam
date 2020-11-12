
% simulation data
% close all;


clear all;
sigma=0.6053;
nu_u=9.82e-6;
nu_l=9.82e-6;
rho_u=1;
rho_l=1000;
lambda = 100e-6;
k_x=2*pi/lambda;
g=9.81*0;
a_0=2E-6;
u_0=0;

% read simulation data

filename = '.\SampleLog';
delimiter = '\t';
startRow = 3;

formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,0.0,'HeaderLines' ,startRow-1, 'ReturnOnError', true);
fclose(fileID);

lines = 1:numel(dataArray{:,14});

time = dataArray{:, 1}(lines);
volume = dataArray{:, 2}(lines);
cx = dataArray{:, 3}(lines);
cy = dataArray{:, 4}(lines);
cz = dataArray{:, 5}(lines);
c2x = dataArray{:, 6}(lines);
c2y = dataArray{:, 7}(lines);
c2z = dataArray{:, 8}(lines);
xPos = dataArray{:, 9}(lines);
xNeg = dataArray{:, 10}(lines);
yPos = dataArray{:, 11}(lines);
yNeg = dataArray{:, 12}(lines);
zPos = dataArray{:, 13}(lines);
zNeg = dataArray{:, 14}(lines);

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

omega0 = sqrt(sigma/(rho_l+rho_u)*power(k_x,3));

tau_offset = -0.0;

t_start=1e-12;
t_end=max(time)+tau_offset/omega0;
N_t=1000;

t=[t_start:(t_end-t_start)/(N_t-1):t_end]';

result=invlap('prosperetti',t,0,1e-20,sigma,nu_u,nu_l,rho_u,rho_l,k_x,g,a_0,u_0);

% figure;
hold on;

plot(time*omega0+tau_offset,(zPos+cz)./lambda,'color','red','linestyle','-');
plot(t*omega0,result/lambda,'color','black','linestyle','-');
