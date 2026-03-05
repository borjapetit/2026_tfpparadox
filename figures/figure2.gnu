reset

set fit quiet

set style line 1 lc rgb '#006400' lw 3
set style line 2 lc rgb 'red'       lw 3 dt 4
set style line 3 lc rgb 'blue'      lw 3 dt 4

hugo0  = "results/h_base_base.txt"
hugo1  = "results/h_einw_base.txt"
hugo2  = "results/h_nok_base.txt"

lucas0 = "results/l_base_base.txt"
lucas2 = "results/l_nok_base.txt"

##########################################################################################

set terminal pdfcairo size 10,5 color font 'arial,20' lw 2

set output 'figures/figure2.pdf'

set xzeroaxis
set yzeroaxis

set xrange [0:20]
set yrange [-5.5:5.5]

set multiplot layout 1,2

####################################################

set title "Hopenhayn's economy"
set xlabel 'Subsidy to small firms (in %)'
set key bottom left
plot hugo0 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 1 title "Baseline", \
     hugo1 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 2 title "Entry cost in labor", \
     hugo2 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 3 title "No capital" #, \

####################################################

set title "Lucas' economy"
set xlabel 'Subsidy to small firms (in %)'
unset key
plot lucas0 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 1 title "Baseline", \
     lucas2 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 3 title "No Capital" 

unset multiplot
unset output

####################################################