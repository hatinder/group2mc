filename2E = '..\MC\cmake-build-debug\EULER001.txt';
delimiterIn = ' ';
headerlinesIn = 1;
D2E = importdata(filename2E, delimiterIn,headerlinesIn);

%histogram(D2E.data(:,2));
hold on
errorbar(log10(D2E.data(:,1)),log10(abs(D2E.data(:,2))),D2E.data(:,3));
plot(log10(D2E.data(:,1)),log10(D2E.data(:,1)));
hold off
xlabel("\Deltat")
ylabel('\textbf{$$log_{10}(|\hat{u}(t,x) - u(t,x)|)$$}','Interpreter','Latex');
title('Weak Euler \Deltat Convergence')
legend('Sim Size: 10^7', 'O(\Deltat)', 'Location', 'bestoutside')
%plot(log10(D2E.data(:,1)),log10(D2E.data(:,1)));
grid on;
grid on;
