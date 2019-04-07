filename2E = 'E:\edubitbucket\group2mc\MC\cmake-build-debug\ST001.txt';
delimiterIn = ' ';
headerlinesIn = 1;
D2E = importdata(filename2E, delimiterIn,headerlinesIn);

histogram(D2E.data(:,1));
%plot(D2E.data(:,1),D2E.data(:,201));
grid on;
