filename2E = '..\MC\cmake-build-debug\RNG002.txt';
delimiterIn = ' ';
headerlinesIn = 1;
D2E = importdata(filename2E, delimiterIn,headerlinesIn);

histogram(D2E.data(:,1));
disp(length(D2E.data(:,1)));
disp(length(D2E.data(:,1))-length(unique(D2E.data(:,1))));
%plot(D2E.data(:,1),D2E.data(:,201));
grid on;
