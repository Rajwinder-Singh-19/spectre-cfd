files = ../spectre_lib/domain.f90 ../spectre_lib/comm.f90 ../spectre_lib/fft.f90 ../spectre_lib/field.f90 ../spectre_lib/turbulence.f90 ../spectre_lib/poisson.f90 ../spectre_lib/navier.f90 ../spectre_lib/diagnostics.f90 ../spectre_lib/vti_writer.f90 ../spectre_lib/solver.f90
modules = domain.mod comm.mod fft.mod field.mod turbulence.mod poisson.mod navier.mod diagnostics.mod vti_writer.mod solver.mod 
fftw_include = -I/usr/include -lfftw3

main: $(files) main.f90
	mkdir output
	mpif90 -o spectre $(files) main.f90 -O3 -ffree-line-length-none $(fftw_include) 
clean: 
	rm -rf output
	rm *.mod
	rm spectre
init: 
	rm -rf output
	mkdir output
