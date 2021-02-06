time = "456.95893"

set term pdf
set output "Ux".time.".pdf"
set logscale x

plot time."/Ux.xy" u 1:2 t "OF" , [0.001:10] x , [10:1000] (1/0.41)*log(x) + 5.2 