lte_dl:
	g++ -c -o lte_dl.o lte_dl.cpp
	g++ -o lte lte_dl.o -lfftw3
