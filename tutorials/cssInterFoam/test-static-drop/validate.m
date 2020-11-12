% close all;

filename = '.\SampleLog';
delimiter = '\t';
startRow = 3;

formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';


fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);

time = dataArray{:, 1};
volume = dataArray{:, 2};
cx1 = dataArray{:, 3};
cy1 = dataArray{:, 4};
cz1 = dataArray{:, 5};
c2x1 = dataArray{:, 6};
c2y1 = dataArray{:, 7};
c2z1 = dataArray{:, 8};
xPos1 = dataArray{:, 9};
xNeg1 = dataArray{:, 10};
yPos1 = dataArray{:, 11};
yNeg1 = dataArray{:, 12};
zPos1 = dataArray{:, 13};
zNeg1 = dataArray{:, 14};
xVel = dataArray{:, 15};
yVel = dataArray{:, 16};
zVel = dataArray{:, 17};
contArea = dataArray{:, 18};
pMeanIn = dataArray{:, 19};
pMeanOut = dataArray{:, 20};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

sigma = 0.4; % surface tension

r0 = power(2*volume*3/pi/4,1/3);

deltaP = 2*sigma*power(r0,-1);

nStart = 30;

hold on;

plot(time(nStart:end),(pMeanIn(nStart:end)-pMeanOut(nStart:end)).*power(deltaP(nStart:end),-1),'xb');



