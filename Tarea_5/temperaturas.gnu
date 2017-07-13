# set terminal cairolatex pdf size 12.5cm,10.5cm standalone colortext \
	header '\newcommand{\hl}[1]{\setlength{\fboxsep}{0.75pt}\colorbox{white}{#1}}'

set terminal pngcairo size 1200,600
set xlabel 'Cantidad de nodos en x'
set ylabel 'Cantidad de nodos en y'
unset key
set palette rgbformulae 33,13,10
set title "Distribucion de temperatura"
set output 'temperatura.png'
set pm3d map interpolate 0,0
set bmargin at screen 0.12
splot 'temperaturas.dat' matrix
