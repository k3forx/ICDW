set ter tikz createstyle
set xrange [-0.00005:0.0008]
set xl "$1/\\Delta \\tau^2$"
set yl "$\\Delta H$"
set format y "$%1.5F$"
set format x "$%1.5F$"
set xtics 0.0002
set ytics
set mxtics
set mytics
set grid
set key right bottom
set xzeroaxis lt -1
set yzeroaxis lt -1

system("rm fit.log")

  f(x) = a*x + b
  fit f(x) 'nmd10to90_hamildiff.dat' u 1:2:3 via a,b

ffac = FIT_STDFIT
system(sprintf("echo 0 %24.15E %24.15E > hamildiff_elp.dat",b,b_err/ffac))

  set output sprintf("hamildiff_elp.tex")
  pl 'nmd10to90_hamildiff.dat' u 1:2:3 w yerrorbars pt 6 ti 'HMC result',\
      'hamildiff_elp.dat' u 1:2:3 w yerrorbars pt 7 ti 'expolation',f(x) w l ti 'fit'
  undefine f,a,b,ffac

system("mv fit.log > elp_fit.log")
