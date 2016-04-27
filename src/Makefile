SHELL = /bin/sh

EXEC = sc_wl

$(shell ./mkFilepath)
$(shell ./mkSrcfiles)
$(shell ./mkDepends Filepath Srcfiles > Depends)

F90 = /usr/bin/gfortran -c
LD = /usr/bin/gfortran
FLAGS := -g -fdefault-real-8
#FLAGS := -g 

INC_NETCDF = /usr/local/include
LIB_NETCDF = /usr/local/lib64
#
#LIB_HDF = ${TACC_HDF5_LIB}


FLAGS += -I${INC_NETCDF}
#LDFLAGS = ${FLAGS} -L${LIB_NETCDF} -L${LIB_HDF} -lnetcdff -lhdf5_hl -lhdf5 -lz
LDFLAGS = ${FLAGS} -L${LIB_NETCDF} -lnetcdff -lnetcdf  

#src_files := $(shell cat Srcfiles | sed s/printnml.f90//g | sed s/make_sounding.F90//g)
src_files := $(shell cat Srcfiles)
obj_files := $(addsuffix .o, $(basename $(src_files)))

.SUFFIXES:
.SUFFIXES: .f .F90 .o

../$(EXEC): $(obj_files)
	$(LD) -o $@ $(obj_files) $(LDFLAGS)

Depends: Srcfiles Filepath
	./mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	./mkSrcfiles

%.o: %.f90
	${F90} ${LDFLAGS} $< -o $@

ifndef nodep
   include Depends
endif


.PHONY: clean
clean:
	rm -f *.mod *.o
