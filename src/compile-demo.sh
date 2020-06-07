g++ -O2 rismhi3d.cpp -o rismhi3d -I$HOME/bin/fftw-3.3.3/include -lm -lz -L$HOME/bin/fftw-3.3.3/lib -lfftw3 -lpthread
g++ -O2 gensolvent.cpp -o gensolvent -I$HOME/bin/fftw-3.3.3/include -lm -lz -L$HOME/bin/fftw-3.3.3/lib -lfftw3 -lpthread
g++ -O2 gmxtop2solute.cpp -o gmxtop2solute -I$HOME/bin/fftw-3.3.3/include -lm -lz -L$HOME/bin/fftw-3.3.3/lib -lfftw3 -lpthread
g++ -O2 ts4sdump.cpp -o ts4sdump -I$HOME/bin/fftw-3.3.3/include -lm -lz -L$HOME/bin/fftw-3.3.3/lib -lfftw3 -lpthread
