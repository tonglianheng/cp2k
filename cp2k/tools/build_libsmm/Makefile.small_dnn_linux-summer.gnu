DIMS         = 1 4 5 6 8 9 13 16 17 22 23 24 26 32 
DIMS_INDICES = $(foreach m,$(DIMS),$(foreach n,$(DIMS),$(foreach k,$(DIMS),$m_$n_$k)))

SI = 1
EI = $(words $(DIMS_INDICES))
INDICES = $(wordlist $(SI),$(EI),$(DIMS_INDICES))

OUTDIR=output_linux-summer.gnu

SRCFILES=$(patsubst %,small_find_%.f90,$(INDICES)) 
EXEFILES=$(patsubst %,$(OUTDIR)/small_find_%.x,$(INDICES)) 
OUTFILES=$(patsubst %,$(OUTDIR)/small_find_%.out,$(INDICES)) 

all: bench 

source: $(SRCFILES) 

%.f90: 
	 .././small_gen.x `echo $* | awk -F_ '{ print $$3" "$$4" "$$5 }'` 1 1 16 ../tiny_gen_optimal_dnn_linux-summer.gnu.out > $@

compile: $(EXEFILES) 

$(OUTDIR)/%.x: %.f90 
	 gfortran -O2 -funroll-loops -ffast-math -ftree-vectorize -march=native -fno-inline-functions $< -o $@ -L/home/tong/opt/acml-5.3.1-gfortran-64bit/gfortran64/lib -Wl,--start-group -lacml -Wl,--end-group 

bench: $(OUTFILES) 

$(OUTDIR)/%.out: $(OUTDIR)/%.x 
	 ./$< > $@

