filename2E = '..\MC\cmake-build-debug\MCCP001.txt';
delimiterIn = ' ';
headerlinesIn = 1;
D2E = importdata(filename2E, delimiterIn,headerlinesIn);

%histogram(D2E.data(:,2));
plot(D2E.data(:,1),D2E.data(:,2));
grid on;
