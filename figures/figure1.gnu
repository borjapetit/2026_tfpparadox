reset

set fit quiet

set style line 1 lc rgb '#006400' lw 3
set style line 2 lc rgb 'red'       lw 3 dt 4
set style line 3 lc rgb 'blue'      lw 3 dt 4

hugo  = "results/h_base_base.txt"
lucas = "results/l_base_base.txt"

##########################################################################################

set terminal pdfcairo size 10,5 enhanced color font 'arial,20' lw 2

set output 'figures/figure1.pdf'
#set output 'figures/figure1.tex'

set xzeroaxis
set yzeroaxis

set xrange [0:20]
set yrange [-5.5:5.5]

set multiplot layout 1,2

####################################################

set title "Hopenhayn's economy"
set xlabel 'Subsidy to small firms (in %)'
set key bottom left
plot hugo every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 1 title "TFP", \
     hugo every ::0::20 using 'atau':'DPhi' smooth bezier with line ls 2  title "Misallocation", \
     hugo every ::0::20 using 'atau':(column("DM") + column("DF")) smooth bezier with line ls 3 title "Selection + Scale"

####################################################

set title "Lucas' economy"
set xlabel 'Subsidy to small firms (in %)'
unset key
plot lucas every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 1, \
     lucas every ::0::20 using 'atau':'DPhi' smooth bezier with line ls 2  , \
     lucas every ::0::20 using 'atau':(column("DM") + column("DF")) smooth bezier with line ls 3

unset multiplot
unset output

####################################################