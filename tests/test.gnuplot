set terminal png white nocrop enhanced size 1024,320 font "arial,8"
set output 'BP_1.png'
set key bmargin center horizontal Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid
set style data lines
set title "Filter BP" 
set xrange [ 0.00000 : 120.00000 ] noreverse nowriteback
x = 0.0
plot 'using.dat' using 1:3 title "power"  

