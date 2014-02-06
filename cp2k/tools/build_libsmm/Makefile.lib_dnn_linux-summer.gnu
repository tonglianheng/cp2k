DIMS         = 1 4 5 6 8 9 13 16 17 22 23 24 26 32 
INDICES = $(foreach m,$(DIMS),$(foreach n,$(DIMS),$(foreach k,$(DIMS),$m_$n_$k)))

OUTDIR=output_linux-summer.gnu

DRIVER=smm_dnn.o 

SRCFILES=$(patsubst %,smm_dnn_%.f90,$(INDICES)) 
OBJFILES=$(patsubst %,$(OUTDIR)/smm_dnn_%.o,$(INDICES)) 

.PHONY: $(OUTDIR)/$(DRIVER) 

all: archive 

source: $(SRCFILES) 

%.f90: 
	 .././lib_gen.x `echo $* | awk -F_ '{ print $$3" "$$4" "$$5 }'` 1 1 16 ../small_gen_optimal_dnn_linux-summer.gnu.out ../tiny_gen_optimal_dnn_linux-summer.gnu.out > $@

compile: $(OBJFILES) 

$(OUTDIR)/%.o: %.f90 
	 gfortran -O2 -funroll-loops -ffast-math -ftree-vectorize -march=native -fno-inline-functions -c $< -o $@ 

$(OUTDIR)/$(DRIVER): 
	 gfortran -O2 -funroll-loops -ffast-math -ftree-vectorize -march=native -fno-inline-functions -c $(notdir $*).f90 -o $@ 

archive: ../lib/libsmm_dnn_linux-summer.gnu.a 

../lib/libsmm_dnn_linux-summer.gnu.a: $(OBJFILES) $(OUTDIR)/$(DRIVER) 
	 ar -r $@ $^ 
	 @echo 'Library produced at /home/tong/tmp/cp2k/cp2k/tools/build_libsmm/lib/libsmm_dnn_linux-summer.gnu.a'

