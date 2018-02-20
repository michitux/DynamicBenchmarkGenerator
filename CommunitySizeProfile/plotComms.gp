set terminal pdf enhanced font 'Helvetica, 16'
set key outside

set xlabel 'Timestep'
set ylabel 'Community Size'
set y2label 'Number of communities'


set y2range [0:]

set ytics nomirror
set y2tics nomirror


set output "commSizeProfile_10K_2_0.9.pdf"
plot "commSizeProfile.log" using 1:2 title 'Min' with lines,\
     	'' using 1:3 title 'Max' with lines, \
     	'' using 1:4 title 'Median' with lines, \
     	'' using 1:5 title 'Avg' with lines, \
     	'' using 1:6 axes x1y2 title '# comms' with lines 

unset output
