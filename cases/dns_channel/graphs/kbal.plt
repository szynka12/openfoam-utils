time ="1599.356255"

set term pdf
set output "kbal".time.".pdf"

set xrange [0:100]

plot \
  time."/P_ii.xy" u 1:($2/2) t "Production, OF" w l, \
  "chan395.kbal" u 2:4 t "chan395" with lines lc rgb "black", \
  \
  time."/E_ii.xy" u 1:(-$2/2) t "Dissipation, OF" w l, \
  "chan395.kbal" u 2:3 notitle with lines lc rgb "black", \
  \
  time."/T_ii.xy" u 1:(-$2/2) t "Diffusion, turbulent, OF" w l, \
  "chan395.kbal" u 2:7 notitle with lines lc rgb "black", \
  \
  time."/D_ii.xy" u 1:($2/2) t "Diffusion, molecular, OF" w l, \
  "chan395.kbal" u 2:8 notitle with lines lc rgb "black", \
  \
  time."/Pi_ii.xy" u 1:($2/2) t "Vel. p. grad. corr., OF" w l, \
  "chan395.kbal" u 2:6 notitle with lines lc rgb "black"
