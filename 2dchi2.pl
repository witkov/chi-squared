# 2dchi2.pl version 1/3/20
# This script does chi-square curve fitting to a 2-parameter linear model
# y = Ax + B
# Three arrays are needed:
    # $x is an array of mean values for the independent variable
    # $y is an array of mean values for the depednent variable
    # $y_err is an array of standard errors (i.e., SD/(sqrt of N)) for y
# Note that this script only handles errors on the dependent (y) variable.
# SOFTWARE DEPENDENCIES: Perl Data Language (PDL) and Gnuplot
# ---------------------------------------------------------------------------
# Copyright (C) 2020 Carey Witkov 
# 
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program. 
# If not, see https://www.gnu.org/licenses/.
#
# -----------------------------------------------------------------------------


use PDL;
use PDL::Ufunc;
use PDL::NiceSlice;
use PDL::Graphics::Prima::Simple;

# Replace sample data with your own data
$x=pdl q[0.10 0.60 1.10 1.60 2.10 2.60 3.10 3.60 4.10 4.60 5.10 5.60 6.10 6.60 7.10 7.60 8.10 8.60 9.10 9.60];  
$y=pdl q[34.1329 98.7892 121.0725 180.3328 236.3683 260.5684 320.9553 380.3028 407.3759 453.7503 506.9685 576.329 602.0845 699.0915 771.2271 796.6707 787.8436 877.0763 915.3649 1000.7312];
$yerr=pdl q[5.8423 9.9393 11.0033 13.4288 15.3743 16.1421 17.9152 19.5014 20.1836 21.3014 22.516 24.0194 24.5374 26.4403 27.771 28.2254 28.0686 29.6155 30.255 31.6343]; 

# calculate sums needed to obtain chi-square
$S1=sum($y**2/$yerr**2);
$S2=sum($x**2/$yerr**2);
$S3=sum(1/$yerr**2);
$S4=sum(($y * $x)/$yerr**2);
$S5=sum($y/$yerr**2);
$S6=sum($x/$yerr**2);

# create parameter grid
$a = zeros(100)->xlinvals(85,110);
$b = zeros(100)->xlinvals(10,35);
$A = zeros(100,100)->xlinvals(85,110);
$B = zeros(100,100)->ylinvals(10,35);

# create empty chi-square grid
$chi2=zeros(100,100);

# calculate chi-square grid over parameters
$chi2=($S1) + ($A**2)*($S2) + ($B**2)*($S3) - 2*$A*$S4 - 2*$B*$S5 + 2*$A*$B*$S6;

# extract chi=square min and indices
$minchi2=min($chi2);
($first,$second)=one2nd($chi2,$chi2->flat->minimum_ind);
$Abest=sclr $a($first);
$Bbest=sclr $b($second);

# display results
print 'Chi-square min = ';
printf "%.3g\n",$minchi2;
print 'Number of data points = ';
print nelem($y); print "\n";
print 'Best-fit parameter A = ';
printf "%.3g\n",$Abest;
print 'Best-fit parameter B = ';
printf "%.3g\n",$Bbest;
print 'A index = ';
print $first . "\n";
print 'B index = ';
print $second . "\n";

# data plot
plot(
-data => ds::Pair($x,$y,plotTypes => [ppair::ErrorBars(y_err=>$yerr),ppair::Squares(filled =>0),],),
-func => ds::Func(sub {$Abest*$_[0]+$Bbest},plotTypes => ppair::Lines,linestyle => lp::ShortDash,colors => cl::LightBlue),
title => 'Data and best-fit line',
x => { label   => 'X' },
y => { label   => 'Y' },
);
$chi2=$chi2-$minchi2;
$r0=transpose($b)->append($chi2); # insert a new 0 column;
$aa=zeros(1)->append($a); # insert a 0 in 1st position of $a pdl;
$r1=($aa)->glue(1,$r0); # insert a new 0 row;
wcols transpose($r1),'contour.dat';

$file = 'file';
$message = <<"END_MSG";
set view map 
unset surface
set contour base
set cntrparam levels discrete 0.001, 2.3, 6.17
set title "1 and 2 sigma contours in A-B plane" font ",14" tc rgb "#606060"
set xlabel "Parameter A"
set ylabel "Parameter B"
splot "contour.dat" nonuniform matrix using 1:2:3 with linespoints pt -1
END_MSG
open $file, '>', 'chisq.plt';
print $file $message;
close $file;
system("chisq.plt");









