#set terminal cairolatex pdf size 12cm,12cm standalone colortext \
#	header '\newcommand{\hl}[1]{\setlength{\fboxsep}{0.75pt}\colorbox{white}{#1}}'

set terminal pngcairo size 600,600
set xlabel 'x [cm]'
set ylabel 'y [cm]'
set palette rgbformulae 33,13,10
set xrange [-2:50]
set yrange [-2:50]
unset key
set title "Calor en [cal/s] para paso de 6.25 [cm]"
set output 'calor2.png'
set bmargin at screen 0.13
plot 'flujo.dat' using 1:2:(5*$3/$5):(5*$4/$5):5 w vector lc palette
