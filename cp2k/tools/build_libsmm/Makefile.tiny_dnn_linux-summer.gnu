DIMS         = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 
DIMS_INDICES = $(foreach m,$(DIMS),$(foreach n,$(DIMS),$(foreach k,$(DIMS),$m_$n_$k)))

SI = 1
EI = $(words $(DIMS_INDICES))
INDICES = $(wordlist $(SI),$(EI),$(DIMS_INDICES))

OUTDIR=output_linux-summer.gnu

SRCFILES=$(patsubst %,tiny_find_%.f90,$(INDICES)) 
EXEFILES=$(patsubst %,$(OUTDIR)/tiny_find_%.x,$(INDICES)) 
OUTFILES=$(patsubst %,$(OUTDIR)/tiny_find_%.out,$(INDICES)) 

all: bench 

source: $(SRCFILES) 

%.f90: 
	 .././tiny_gen.x `echo $* | awk -F_ '{ print $$3" "$$4" "$$5 }'` 1 1 > $@

compile: $(EXEFILES) 

$(OUTDIR)/%.x: %.f90 
	 gfortran -O2 -funroll-loops -ffast-math -ftree-vectorize -march=native -fno-inline-functions $< -o $@  

bench: $(OUTFILES) 

$(OUTDIR)/%.out: $(OUTDIR)/%.x 
	  ./$< > $@ 

