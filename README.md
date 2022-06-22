## LTE_PHY

### The project has implementations of the various downlink channels for 4G LTE in C++.
### The following channels are implemented:

- Physical Broadcast Channel (PBCH)
- Physical Control Format Indicator Channel (PCFICH)
- Physical Hybrid ARQ Indicator Channel (PHICH)
- Physical Downlink Control Channel (PDCCH)

### The various encoding and mappings used in the different channels are also implemented.

### Refer to *lte_channels.h* for the various function defines.

### The other header files have structures defined to hold data from various channels

### The file lte_channels.cpp was run using g++ compiler

### Before running install fftw3 library: sudo apt-get install libfftw3-dev

### Run make lte_channels

### There is an int main() inside lte_channels.cpp which is preconfigured with values.
### This can be changed depending on the bandwidth,Number of antennas,PHICH_Ng_value and cell_id.

### *Document referred for implementation is 3GPP TS 36.211 v10.1.0*
