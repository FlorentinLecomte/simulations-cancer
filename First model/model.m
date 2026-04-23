x = linspace(0,100,1);
t = linspace(0,10,20);
m = 0;
sol = pdepe(m,@modelpde,@modelic,@modelbc,x,t);
n = sol(:,:,1);
f = sol(:,:,2);
m = sol(:,:,3);
surf(x,t,n)
title('n(x,t): Concentration of cancer cells')
xlabel('Distance x')
ylabel('Time t')
saveas(gcf, 'surfn.png');