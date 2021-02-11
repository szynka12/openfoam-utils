time ="1599.356255"

set term pdf
set output "kbal".time.".pdf"
set logscale x

plot time."/P_ii.xy" u 1:2 t "Production, OF" ,\
  "chan395.kbal" u 2:4 t "Production, chan395"
