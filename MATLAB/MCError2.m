filename2E = '..\MC\cmake-build-debug\MCERROR001.txt';
delimiterIn = ' ';
headerlinesIn = 1;
D2E = importdata(filename2E, delimiterIn,headerlinesIn);

%histogram(D2E.data(:,2));
hold on
errorbar(log10(D2E.data(:,1)),log10(D2E.data(:,2)),D2E.data(:,3));
plot(log10(D2E.data(:,1)),log10(D2E.data(:,1)));
hold off
grid on;
