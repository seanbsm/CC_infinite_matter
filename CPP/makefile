# the compiler to use
CPP  := g++
FORTC := gfortran 

# -no-pie -g -pg for gprof

TARGET := run

# Source, object, and dependency files for c++, cuda, and fortran
SRC_FILES := $(shell find ./ -regex [^\#]*\\.cpp$)
OBJ_FILES := $(SRC_FILES:.cpp=.o)
DEP_FILES := $(SRC_FILES:.cpp=.d)

SRCFORT_FILES := $(shell find ./ -regex [^\#]*\\.f90$)
OBJFORT_FILES := fort_run.o
OBJFORT_FILES := $(SRCFORT_FILES:.f90=.o)
#~ OBJFORT_FILES := $(SRCFORT_FILES:./Potentials/chp/%.f90=%.o)

# -ggdb
# Flags
CPPFLAGS := -Wall -std=c++17 -O2 -DMKL_ILP64 -m64 -DMKL_Complex16="std::complex<double>" -DMKL_Complex8="std::complex<float>" -I$(MKLROOT)/include
FORTFLAGS := -O2 -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC

# Libraries
LDLIBS := -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -lgfortran -lgsl -lgslcblas


# Symbol explenations
# @ = recipe name (so $(TARGET) in this case)
# ^ = all dependencies (so ALL $(OBJ_FILES))
# < = first dependency (so first of OBJ_FILES)

# Include the .d files
-include $(DEP_FILES)
#~ printfortfiles:
#~ 	echo $(SRCFORT_FILES)
#~ 	echo $(OBJFORT_FILES)

# Make executable
$(TARGET): $(OBJ_FILES) $(OBJFORT_FILES)
	$(CPP) $^ -o $@ $(LDFLAGS) $(LDLIBS)

# Make object files
%.o: %.cpp makefile
	$(CPP) $(CPPFLAGS) -MMD -MP -c $< -o $@
	
$(OBJFORT_FILES): $(SRCFORT_FILES)
	$(FORTC) $(FORTFLAGS) -c $^

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJ_FILES) $(OBJFORT_FILES) $(DEP_FILES) *.o *.mod *.exe

.PRECIOUS: $(SRCFORT_FILES)
