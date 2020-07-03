# 1dchi2.pl version 1/3/20
# This script does chi-square curve fitting to a 1-parameter linear model
# y = Ax
# Needed are three arrays:
    # $x is an array of mean values for the independent variable
    # $y is an array of mean values for the depednent variable
    # $y_err is an array of standard errors (i.e., SD/(sqrt of N)) for y
# Note that this script only handles errors on the dependent (y) variable.
# SOFTWARE DEPENDENCIES: Perl Data Language (PDL)
# ---------------------------------------------------------------------------
# AUTHOR:  Carey Witkov, witkov@fas.harvard.edu, and Keith Zengel, zengel@fas.harvard.edu
# COPYRIGHT: This script is copyright (2017) Carey Witkov and Keith Zengel.
# LICENSE:  Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
# https://creativecommons.org/licenses/by-nc/4.0/
# DISCLAIMER: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ---------------------------------------------------------------------------

use PDL;
use PDL::Graphics::Prima::Simple;
use PDL::Ufunc;
use PDL::NiceSlice;

# Replace sample data with your own data
$x=pdl q[0.00873 0.01765 0.02619 0.03502 0.04375];  
$y=pdl q[1.13806 2.35315 3.93546 4.80486 5.50653];
$yerr=pdl q[0.0659 0.0876 0.1080 0.1121	0.2025]; 

# calculate sums needed to obtain chi-square
$S1=sum($y**2/$yerr**2);
$S2=sum($x**2/$yerr**2);
$S4=sum(($y * $x)/$yerr**2);

# create parameter grid
$a = zeroes(100)->xlinvals(133,141);

# calculate chi-square over parameter grid
$chi2 = $S1 + ($a**2)*($S2) - 2*$a*$S4;  # calculate chi-square using sum arrays
     
# extract chi-square min and indices
$min_val = minimum($chi2);
$idx = minimum_ind($chi2);
$minchi2=$min_val;

# extract best parameter values at indices of chi-square min
$Abest=$a($idx);

# display results
print 'Chi-square min = ';
printf "%.3g\n",$minchi2;
print 'Number of data points = ';
print nelem($y); print "\n";
print 'Best-fit parameter A = ';
$Abest = sclr($Abest);
printf "%.3g\n",$Abest;

# plot of data and model line
plot(
-data => ds::Pair($x,$y,plotTypes => [ppair::ErrorBars(y_err=>$yerr),ppair::Squares(filled =>0),],),
-func => ds::Func(sub {$Abest*$_[0]+$Bbest},plotTypes => ppair::Lines,linestyle => lp::ShortDash,colors => cl::LightBlue),
title => 'Data and best-fit slope',
x => { label   => 'X' },
y => { label   => 'Y' },
);

# plot of chi-square vs. parameter A (slope)
plot(
# the experimental data
-data => ds::Pair($a,$chi2,plotTypes => [ppair::Blobs(radius=>2,colors => cl::LightBlue)],),
title => 'Chi-square vs. A (slope)',
x => { label   => 'Parameter A (slope)' },
y => { label   => 'Chi-square' },
);


