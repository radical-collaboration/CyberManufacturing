clc
clear
close all

d = importdata('d50OverTime.csv');
d = d.data();
d(2,5) = 0.0;
figure()
plot(d(:, 2), d(:, 3),'r')
hold on
plot(d(:, 2), d(:, 4),'c')
hold on
plot(d(:, 2), d(:, 5),'b')
hold on
plot(d(:, 2), d(:, 6),'k')
legend('comp1', 'comp2', 'comp3', 'comp4')
title('d50')
hold off

d = importdata('d10OverTime.csv');
d = d.data();
d(2,5) = 0.0;
figure()
plot(d(:, 2), d(:, 3),'r')
hold on
plot(d(:, 2), d(:, 4),'c')
hold on
plot(d(:, 2), d(:, 5),'b')
hold on
plot(d(:, 2), d(:, 6),'k')
legend('comp1', 'comp2', 'comp3', 'comp4')
title('d10')
hold off

d = importdata('d90OverTime.csv');
d = d.data();
d(2,5) = 0.0;
figure()
plot(d(:, 2), d(:, 3),'r')
hold on
plot(d(:, 2), d(:, 4),'c')
hold on
plot(d(:, 2), d(:, 5),'b')
hold on
plot(d(:, 2), d(:, 6),'k')
legend('comp1', 'comp2', 'comp3', 'comp4')
title('d90')
hold off
