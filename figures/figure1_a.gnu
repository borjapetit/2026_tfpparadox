reset

set fit quiet

set style line 1 lc rgb '#006400' lw 3
set style line 2 lc rgb 'red'       lw 3 dt 4
set style line 3 lc rgb 'blue'      lw 3 dt 4

hugo  = "results/h_base_base.txt"
lucas = "results/l_base_base.txt"

##########################################################################################

set terminal pdfcairo size 10,14 enhanced color font 'arial,20' lw 2

set output 'figures/figure1_a.pdf'

set xzeroaxis
set xrange [0:21]

set multiplot layout 3,2

set title "Output"
set xlabel 'Subsidy to small firms (in %)'
set key bottom left
plot hugo  every ::1::20 using 'atau':'aY' smooth bezier with line ls 1 dt 1 title "Hopenhayn", \
     lucas every ::1::20 using 'atau':'aY' smooth bezier with line ls 2 dt 1 title "Lucas"

set title "Consumption"
set xlabel 'Subsidy to small firms (in %)'
unset key
plot hugo  every ::1::20 using 'atau':'aC' smooth bezier with line ls 1 dt 1 title "Hopenhayn", \
     lucas every ::1::20 using 'atau':'aC' smooth bezier with line ls 2 dt 1 title "Lucas"

set title "Wages"
set xlabel 'Subsidy to small firms (in %)'
unset key
plot hugo  every ::1::20 using 'atau':'wage' smooth bezier with line ls 1 dt 1 title "Hopenhayn", \
     lucas every ::1::20 using 'atau':'wage' smooth bezier with line ls 2 dt 1 title "Lucas"

set title "Labor"
set xlabel 'Subsidy to small firms (in %)'
unset key
plot hugo  every ::1::20 using 'atau':'aN' smooth bezier with line ls 1 dt 1 title "Hopenhayn", \
     lucas every ::1::20 using 'atau':'aN' smooth bezier with line ls 2 dt 1 title "Lucas"

set title "Mass of firms"
set xlabel 'Subsidy to small firms (in %)'
unset key
plot hugo  every ::1::20 using 'atau':'M' smooth bezier with line ls 1 dt 1 title "Hopenhayn", \
     lucas every ::1::20 using 'atau':'M' smooth bezier with line ls 2 dt 1 title "Lucas"

set title "Average productivity"
set xlabel 'Subsidy to small firms (in %)'
unset key
plot hugo  every ::1::20 using 'atau':'Etau0' smooth bezier with line ls 1 dt 1 title "Hopenhayn", \
     lucas every ::1::20 using 'atau':(-column("atau")) smooth bezier with line ls 2 dt 1 title "Lucas"


unset multiplot
unset output

####################################################