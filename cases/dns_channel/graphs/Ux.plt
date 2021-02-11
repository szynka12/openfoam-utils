time ="1599.356255"

set term pdf
set output "Ux".time.".pdf"
set logscale x

plot time."/Ux.xy" u 1:2 t "OF" ,\
  [0.001:10] x notitle w l lc rgb "red",\
  [10:1000] (1/0.41)*log(x) + 5.2 notitle w l lc rgb "red" 
