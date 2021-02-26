time ="913.91786286248521"
target_time = "../constant/targetMesh/graphs/".time
target_time_uni = "../constant/targetMesh_uniform/graphs/".time

set term pdf
set output "Ux".time.".pdf"
set logscale x

plot time."/Ux.xy" u 1:2 t "OF" w l lc rgb "black" ,\
  target_time."/Ux.xy" u 1:2 t "OF <d. avg>" w p pointtype 7 ps 0.4 ,\
  target_time_uni."/Ux.xy" u 1:2 t "OF <d. avg>" w p pointtype 7 ps 0.4 ,\
  [0.001:10] x notitle w l lc rgb "red",\
  [10:1000] (1/0.41)*log(x) + 5.2 notitle w l lc rgb "red" 
  
