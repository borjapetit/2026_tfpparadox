
# MAKEFILE for the recplication code

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
pwd = /mnt/c/Users/bpetit/Dropbox/research/projects/2021_bunching/code_tfp
name = model
endif

ifeq ($(UNAME_S),Darwin)
pwd = /Users/borjapetit/Library/CloudStorage/Dropbox/research/projects/2021_bunching/code_tfp
name = model
endif

fcomp = gfortran
flags = -fopenmp -finit-local-zero -ffixed-line-length-450 -J $(pwd)/compiledfiles/
files = toolkit.f90 parameters.f90 lucas.f90 hugo.f90 main.f90

model: $(files)
	rm -f ${name}
	rm -f model_check
	rm -f */*.mod
	rm -f */*.mod0
	rm -f *.mod
	rm -f *.mod0
	rm -f *.log
	rm -rf *.dSYM
	echo "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
	$(fcomp)  -O3 ${flags} ${files} -o ${name}

check: $(files)
	rm -f model*
	rm -f */*.mod
	rm -f */*.mod0
	rm -f *.mod
	rm -f *.mod0
	rm -f *.logm
	echo "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
	$(fcomp) -fcheck=all -fbacktrace -Wall -g -O2 ${flags} ${files} -o test

clean:
	rm -f model*
	rm -f */*.mod
	rm -f */*.mod0
	rm -f *.mod
	rm -f *.mod0
	rm -f *.log


# End of the makefile

