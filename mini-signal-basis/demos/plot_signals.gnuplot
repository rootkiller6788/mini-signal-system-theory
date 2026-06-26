# Gnuplot script to visualize elementary signals
set terminal png size 800,600
set output 'signals.png'
set title 'Elementary Signals'
set xlabel 'Time (s)'
set ylabel 'Amplitude'
set grid
plot sin(2*pi*x) title 'Sinusoid',      exp(-3*x) title 'Exponential Decay',      (abs(x)<0.5?1:0) title 'Rect Pulse'
