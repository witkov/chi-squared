clear all
close all

% Measure each set of filters 3 times and record the values
l1 = .72;
l2 = 1.00;
l3 = 1.27;
l4 = 1.5;
l5 = 1.65;
density = 0.0873;
length = [l1 l2 l3 l4 l5]
m = density*length;
g=9.81;

%measure the terminal velocity for each set of filters 5 times
%and record the values
F1 = [1.9859 1.5292 1.7136 2.4857 2.3128];
F2 = [3.9875 3.1660 2.7310 3.0544 2.6344];
F3 = [3.1135 4.4017 3.6917 3.6077 3.1158];
F4 = [4.1339 4.0360 4.2656 3.8470 3.4142];
F5 = [4.5564 3.4852 4.2987 4.7974 4.0171];

F = [mean(F1) mean(F2) mean(F3) mean(F4) mean(F5)];
sigma_v=[std(F1) std(F2) std(F3) std(F4) std(F5)]/sqrt(5)

% x,y data arrays and y-error array
x = m*g;
y = F; 
yerr = sigma_v; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% DO NOT EDIT BELOW THIS LINE %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate sums needed to obtain chi-square
s_yy=sum(y.^2./yerr.^2);
s_xx=sum(x.^2./yerr.^2);
s_xy=sum((y.*x)./yerr.^2);
s_y =sum(y./yerr.^2);
s_x =sum(x./yerr.^2);
s_0 =sum(1./yerr.^2);

% by completing the square, we rewrite chi-squared as
% sum((y_i - A x_i - B)^2/sigma_i^2 
% = (A-A^*)^2/sigma_A^2
% + (B-B^*)^2/sigma_B^2 
% + 2*rho*(A-A^*)(B-B^*)/sigma_A*Sigma_A
% + \chi^2_{min}
A_best = (s_0*s_xy - s_x*s_y)/(s_0*s_xx - s_x^2)
sigma_A = 1/sqrt(s_xx);
B_best = (s_y*s_xx - s_x*s_xy)/(s_0*s_xx - s_x^2)
sigma_B = 1/sqrt(s_0);
rho = s_x/sqrt(s_xx*s_0);
minchi2 = (s_0*s_xy^2 - 2*s_x*s_y*s_xy + s_y^2*s_xx)/(s_x^2 - s_0*s_xx) + s_yy

%define chi-squared plot interval
A_interval = sqrt(6.17*sigma_A^2/(1-rho^2));
B_interval = sqrt(6.17*sigma_B^2/(1-rho^2));
AB_interval = [B_best-1.1*B_interval B_best+1.1*B_interval A_best-1.1*A_interval A_best+1.1*A_interval];

% Contour plot
figure(1);
paraboloid = @(B,A) (A-A_best).^2/sigma_A.^2 + (B-B_best).^2./sigma_B.^2  + 2*rho*(A-A_best).*(B-B_best)./(sigma_A.*sigma_B);
fcontour(paraboloid, AB_interval,'LevelList',[2.3 6.17],'MeshDensity',200);
xlabel('B (intercept)','Fontsize',16);
ylabel('A (slope)','Fontsize',16);
title('68% and 95% contours of \chi^2 in the AB plane','Fontsize',16);
hold on;
plot(B_best,A_best, 'x');
hold off;

%define fit plot interval around x_min and x_max
delta_x  = range(x);
xmintoxmax = [min(x)-0.1*delta_x max(x)+0.1*delta_x];

%plot data and line of best fit
figure(2);
errorbar(x, y, yerr,'x');
xlabel('x','Fontsize',16);
ylabel('y','Fontsize',16);
title('y vs x data and line-of-best-fit','Fontsize',16);
hold on;
fplot(@(x) A_best.*x + B_best, xmintoxmax);
hold off;

n_best = 1/A_best
k_best = exp(-B_best/A_best)


