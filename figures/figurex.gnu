reset

set fit quiet

set style line 1 lc rgb '#006400' lw 3
set style line 2 lc rgb 'red'       lw 3 dt 4
set style line 3 lc rgb 'blue'      lw 3 dt 4

lucas0 = "results/l_base_s_base.txt"
lucas1 = "results/l_base_s_gamma_0.txt"
lucas2 = "results/l_base_s_gamma_1.txt"
lucas3 = "results/l_base_s_xi_0.txt"
lucas4 = "results/l_base_s_xi_1.txt"

##########################################################################################

set terminal epslatex size 10,5 color lw 2 #font 'arial,20'

set output 'figures/figurex.tex'

set xzeroaxis
set yzeroaxis

set xrange [0:20]
set yrange [-5.5:5.5]

set multiplot layout 1,2

####################################################

set title '$\gamma$'
set xlabel 'Subsidy to small firms (in \%)'
set key bottom left
plot lucas0 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 1 title "Baseline", \
     lucas1 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 2 title "Lower", \
     lucas2 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 3 title "Higher" #, \

####################################################

set title "Lucas' economy"
set xlabel 'Subsidy to small firms (in \%)'
unset key
plot lucas0 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 1 title "Baseline", \
     lucas3 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 2 title "Lower" , \
     lucas4 every ::0::20 using 'atau':'DTFP' smooth bezier with line ls 3 title "Higher" 

unset multiplot
unset output

####################################################