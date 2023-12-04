roqj.o: roqj.h roqj.cpp
	g++ roqj.cpp -c -o roqj.o -std=c++20 -O3 -ffast-math -fno-math-errno

enm: Examples/enm.cpp roqj.o
	g++ Examples/enm.cpp roqj.o -o Examples/enm.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/enm.x
	python3.10 Examples/plot.py "Eternally non-Markovian" "tr$$ [\rho\sigma_z] $$ "

enm_x: Examples/enm_x.cpp roqj.o
	g++ Examples/enm_x.cpp roqj.o -o Examples/enm_x.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/enm_x.x
	python3.10 Examples/plot.py "Eternally non-Markovian, jumps to $$ |\pm> $$ " "tr$$ [\rho\sigma_z] $$ "

nM_gamma_z: Examples/nM_gamma_z.cpp roqj.o
	g++ Examples/nM_gamma_z.cpp roqj.o -o Examples/nM_gamma_z.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/nM_gamma_z.x
	python3.10 Examples/plot.py "Non-P-divisible $$ \gamma_z < \sqrt{\gamma_+\gamma_-} $$ " "tr$$ [\rho\sigma_z] $$ "

nM_gamma_pm: Examples/nM_gamma_pm.cpp roqj.o
	g++ Examples/nM_gamma_pm.cpp roqj.o -o Examples/nM_gamma_pm.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/nM_gamma_pm.x
	python3.10 Examples/plot.py "Non-P-divisible $$ \gamma_\pm < 0 $$ " "tr$$ [\rho\sigma_z] $$ "

driven: Examples/driven.cpp roqj.o
	g++ Examples/driven.cpp roqj.o -o Examples/driven.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven.x
	python3.10 Examples/plot.py "Ph cov in x direction, $$ \gamma_\pm = 1 $$ , $$ H = 10 \sigma_z $$ " "tr $$ [\rho\sigma_x] $$ "

driven_enm: Examples/driven_enm.cpp roqj.o
	g++ Examples/driven_enm.cpp roqj.o -o Examples/driven_enm.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_enm.x
	python3.10 Examples/plot.py "ENM x direction, $$ H = 10 \sigma_z $$ " "tr $$ [\rho\sigma_x] $$ "