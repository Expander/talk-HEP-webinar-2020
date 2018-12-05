set term pdf size 4in,4in
set key box bottom right height 0.5
set grid
set logscale x
set format x '10^{%L}'
set xlabel 'M_S / GeV'
set ylabel 'M_h / GeV'
set output 'SplitMSSM_non_degenerate.pdf'

data2  = 'SplitMSSM_degenerate_MS_TB-2_Xt-2.44949_Mlow-1500.dat'
data10 = 'SplitMSSM_degenerate_MS_TB-10_Xt-2.44949_Mlow-1500.dat'
data20 = 'SplitMSSM_degenerate_MS_TB-20_Xt-2.44949_Mlow-1500.dat'
data50 = 'SplitMSSM_degenerate_MS_TB-50_Xt-2.44949_Mlow-1500.dat'

set style line 1 lc rgb '#FF0000' lt 1 lw 0 dt 1
set style line 2 lc rgb '#0000FF' lt 1 lw 0 dt 1

plot \
     data2  u 1:4:5 t 'tan {/Symbol b} = 2'  w filledcurves ls 1 fs transparent solid 0.3, \
     data50 u 1:4:5 t 'tan {/Symbol b} = 50' w filledcurves ls 2 fs transparent solid 0.3
