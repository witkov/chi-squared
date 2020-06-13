% chi2_1D.m version 6/13/20
% This script does chi-square curve fitting to the 1-parameter linear model y = Ax
%
% Three arrays are needed:
%  x is an array of mean values for the independent variable
%  y is an array of mean values for the dependent variable
%  y_err is an array of standard errors (SD/(sqrt of N)) for y
% 
% Note this script only handles errors on dependent (y) variable.
% 
% This script reproduces Figs. 21 and 22 in our book:
% "Chi-Squared Data Analysis and Model Testing for Beginners"
% Oxford University Press, 2019.
% 
%  ---------------------------------------------------------------------------
% Copyright (C) 2020 Carey Witkov and Keith Zengel
% 
% This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation; either version 3 of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. 
% If not, see https://www.gnu.org/licenses/.
%
% -------------------------------------------------------------------------

clear all; 
close all;

% Measure the mass of each set of filters
m1 = 0.89;
m2 = 1.80;
m3 = 2.67;
m4 = 3.57;
m5 = 4.46;
m = [m1 m2 m3 m4 m5]*0.001;
g=9.81;

%measure the terminal velocity for each set of filters 5 times and record the values
v1 = [1.150 1.026 1.014 1.134 1.010];
v2 = [1.443 1.492 1.583 1.563 1.589];
v3 = [2.031 1.939 1.941 2.067 1.941];
v4 = [2.181 2.202 2.199 2.109 2.269];
v5 = [2.507 2.292 2.303 2.364 2.267];

v = [mean(v1) mean(v2) mean(v3) mean(v4) mean(v5)];
sigma_v=[std(v1) std(v2) std(v3) std(v4) std(v5)]/sqrt(5);

% x,y data arrays and y-error array
x = m;
y = v.^2; 
yerr = 2*v.*sigma_v; 

%----------------------------
% DO NOT EDIT BELOW THIS LINE 
%----------------------------

% calculate sums needed to obtain chi-square
s_yy=sum(y.^2./yerr.^2);
s_xx=sum(x.^2./yerr.^2);
s_xy=sum((y.*x)./yerr.^2);

% by completing the square, we rewrite chi-squared as
% sum((y_i - A x_i)^2/sigma_i^2 = (A-A^*)^2/sigma_A^2 + \chi^2_{min}

A_best = s_xy/s_xx
sigma_A = 1/sqrt(s_xx)
minchi2 = s_yy - s_xy^2/s_xx

%define fit plot interval around x_min and x_max
delta_x  = range(x);
xmintoxmax = [min(x)-0.1*delta_x max(x)+0.1*delta_x];

%plot data and line of best fit
figure()
errorbar(x, y, yerr,'x')
xlabel('m [kg]','Fontsize',16)
ylabel('v^2 [m^2/s^2]','Fontsize',16)
hold on
fplot(@(x) A_best.*x, xmintoxmax)
legend('Data','Best Fit','Location','northwest')
hold off

%define chi-squared plot interval as A^* +- 2*sigma_A
twosigma = [A_best-2*sigma_A A_best+2*sigma_A];

%plot chi-squared parabola and minchi+1
figure()
fplot( @(A) (A - A_best).^2/sigma_A.^2 + minchi2, twosigma)
xlabel('A','Fontsize',16)
ylabel('\chi^2','Fontsize',16)
hold on;
fplot(minchi2+1, twosigma)
legend('\chi^2','\chi^2_{min}+1','Location','north')

