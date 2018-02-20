set terminal pdf
set logscale x
set logscale y
set output 'DegreeDynamic.pdf'
set multiplot layout 2,2
set tmargin 1
plot "degreeDist_0" using 1:2 title "T=0" with linespoints linetype 1
#
plot "degreeDist_250" using 1:2 title "T=250" with linespoints linetype 2
#
plot "degreeDist_500" using 1:2 title "T=500" with linespoints linetype 3
#
plot "degreeDist_750" using 1:2 title "T=750" with linespoints linetype 4
#
unset multiplot
#####################
set output 'MembershipDynamic.pdf'
set multiplot layout 2,2
set tmargin 1
plot "commMemDist_0" using 1:2 title "T=0" with linespoints linetype 1
#
plot "commMemDist_250" using 1:2 title "T=250" with linespoints linetype 2
#
plot "commMemDist_500" using 1:2 title "T=500" with linespoints linetype 3
#
plot "commMemDist_750" using 1:2 title "T=750" with linespoints linetype 4
#
unset multiplot
#####################

set output 'CommSizeDynamic.pdf'
set multiplot layout 2,2
set xrange [10:1000]
set tmargin 1
plot "commSizesDist_0" using 1:2 title "T=0" with linespoints linetype 1
#
plot "commSizesDist_250" using 1:2 title "T=250" with linespoints linetype 2
#
plot "commSizesDist_500" using 1:2 title "T=500" with linespoints linetype 3
#
plot "commSizesDist_750" using 1:2 title "T=750" with linespoints linetype 4
#
unset multiplot

