reset
set terminal lua #fulldoc \
#solid originreset plotsize 9,6.5
set output 'plotExecTimeMixedVsPlainBig.tex'

set format x "%g"
set format y "%g"

#unset key
#set label 1 '' #at 5.5, 0.0193014, 0 left norotate back nopoint offset character 0, 0, 0
#set label 2 '$\sigma$' at 7.73607, 0.108213, 0 left norotate back nopoint offset character 0, 0, 0

set xlabel 'Monte Carlo cycles'  
set xrange [40000 : 2050000] noreverse nowriteback
set mxtics 5

set ylabel 'Execution time, s'
set yrange [0 : 30] noreverse nowriteback
set mytics 5

set key 700000,28#top left

plot 'execTimeMixedHeOpt1.data' using 1:2 lt 1 lw 3 w lp t 'Mixed -O1', 'performanceCppOpt1.data' using 2:1 lt 2 lw 3 w lp t 'C++ -O1','execTimeMixedHeOpt2.data' using 1:2 lt 3 lw 3 w lp t 'Mixed -O2','performanceCppOpt2.data' using 2:1 lt 4 lw 5 w lp t 'C++ -O2','execTimeMixedHeOpt3.data' using 1:2 lt 5 lw 3 w lp t 'Mixed -O3','performanceCppOpt3.data' using 2:1 lt 6 lw 3 w lp t 'C++ -O3'












        
