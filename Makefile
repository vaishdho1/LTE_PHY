lte_channels:
	g++ -c -o lte_channels.o lte_channels.cpp
	g++ -o lte_channels lte_channels.o -lfftw3
