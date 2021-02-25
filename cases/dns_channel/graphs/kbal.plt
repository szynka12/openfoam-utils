time ="1599.356255"
time_tM ="1599.356255_targetMesh"


set term pdf
set output "kbal".time.".pdf"

set xrange [0:5]

plot \
  time."/P_ii.xy" u 1:($2/2) t "Production, OF" w l, \
  "chan395.kbal" u 2:4 t "chan395" with lines lc rgb "black" dashtype 7, \
  \
  time."/E_ii.xy" u 1:(-$2/2) t "Dissipation, OF" w l, \
  "chan395.kbal" u 2:3 notitle with lines lc rgb "black" dashtype 7, \
  \
  time."/T_ii.xy" u 1:(-$2/2) t "Diffusion, turbulent, OF" w l, \
  "chan395.kbal" u 2:7 notitle with lines lc rgb "black" dashtype 7, \
  \
  time."/D_ii.xy" u 1:($2/2) t "Diffusion, molecular, OF" w l, \
  "chan395.kbal" u 2:8 notitle with lines lc rgb "black" dashtype 7, \
  \
  time."/Pi_ii.xy" u 1:($2/2) t "Vel. p. grad. corr., OF" w l lc rgb "green", \
  "chan395.kbal" u 2:6 notitle with lines lc rgb "black" dashtype 7 
