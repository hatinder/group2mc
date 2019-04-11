filename2MCE1 = '..\MC\cmake-build-debug\MCAONERR000.txt';
filename2MCE2 = '..\MC\cmake-build-debug\MCAONERR001.txt';
filename2MCE3 = '..\MC\cmake-build-debug\MCAONERR002.txt';
filename2MCE4 = '..\MC\cmake-build-debug\MCAONERR003.txt';
delimiterIn = ' ';
headerlinesIn = 1;
D2MCE1 = importdata(filename2MCE1, delimiterIn,headerlinesIn);
D2MCE2 = importdata(filename2MCE2, delimiterIn,headerlinesIn);
D2MCE3 = importdata(filename2MCE3, delimiterIn,headerlinesIn);
D2MCE4 = importdata(filename2MCE4, delimiterIn,headerlinesIn);

hold on;
histfit(D2MCE1.data(:,1));
%histfit(D2MCE2.data(:,1));
%histfit(D2MCE3.data(:,1));
histfit(D2MCE4.data(:,1));
hold off;

xlabel("Mean")
ylabel("Density")
title('Asset-Or-Nothing Call - Monte Carlo Error Convergence')
legend('Sim Size: 10^3','histfit', 'Sim Size: 4x10^3', 'Location', 'east')
%plot(log10(D2E.data(:,1)),log10(D2E.data(:,1)));
grid on;
