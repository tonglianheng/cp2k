# -*- coding: utf-8 -*-

from math import ceil

class Kernel_dnt_small(object):
    def __init__(self, **params):
        self.__dict__.update(params)
        self.name  = "cusmm_dnt_small_"
        self.name += "_".join([str(params[k]) for k in sorted(params.keys())])
        assert(self.threads * self.minblocks <= 2048)
        min_threads = ((self.m+self.tile_m-1)/self.tile_m) * ((self.n+self.tile_n-1)/self.tile_n)
        assert(min_threads <= self.threads)

    def __repr__(self):
        return("<%s>"%self.name)

    def can_handle(self, m, n, k):
        return(self.m==m and self.n==n and self.k==k)

    def include(self):
        return("cusmm_dnt_small.h")

    def launcher_code(self):
       output  = "int launch_"+self.name+"(int *param_stack, int stack_size, "
       output += "cudaStream_t stream, int m_max, int n_max, int k_max, "
       output += "double *a_data, double *b_data, double *c_data){\n"
       output += "int shared_size = 0;\n"
       output += "//%s\n"%str(self.__dict__)
       output += "int careful = (stack_size / %(grouping)d);\n"%self.__dict__
       output += "int nruns = stack_size - careful * %(grouping)d;\n"%self.__dict__
       output += "typedef void (*kernel)(const int*, int, int, double*, double*, double*);\n"
       output += "static kernel kern_func = cusmm_dnt_small<%(m)d,%(n)d,%(k)d,%(tile_m)d,%(tile_n)d,%(threads)d,%(grouping)d,%(minblocks)d>;\n"%self.__dict__
       output += "static bool configured = false;\n"
       output += "if(configured == false){\n"
       output += "  cudaError_t err = cudaFuncSetSharedMemConfig(kern_func, cudaSharedMemBankSizeEightByte);\n"
       output += "  if(err != cudaSuccess) return(-1);\n"
       output += "  configured = true;\n"
       output += "}\n"
       output += "kern_func<<< ((stack_size + %(grouping)d - 1) / %(grouping)d), %(threads)d, shared_size, stream >>>\n"%self.__dict__
       output += "(param_stack, careful, nruns, \n"
       output += "a_data, b_data, c_data);\n"
       output += "return(0);\n"
       output += "}\n"
       return(output)

    @staticmethod
    def promising_parameters(m, n, k):
        params = []
        grouping = 16
        for minblocks in (1,4,8,12):
            for threads in (64, 96, 128):
                if(threads * minblocks > 2048):
                    continue
                for tm in (1, 2,):
                    for tn in (1, 2,):
                        min_threads = ((m+tm-1)/tm) * ((n+tn-1)/tn)
                        if(min_threads > threads):
                            continue # not enough threads to cover result matrix

                        if(threads > 4*min_threads):
                            continue #heuristic: too many threads unused during calculation

                        if(tm*tn*threads*minblocks > 10000):
                            continue #heuristic: too many registers used

                        buf_sz = max(m*k + k*n, m*n)
                        sizeof_int = 4; sizeof_double = 8
                        smem_tot = buf_sz*sizeof_double + 4*grouping*sizeof_int
                        if(smem_tot*minblocks > 48*1024):
                            continue # uses too much shared memory

                        params.append({'m':m, 'n':n, 'k':k,
                                       'tile_m':tm, 'tile_n':tn,
                                       'threads':threads,
                                       'grouping':grouping,
                                       'minblocks':minblocks})
        return(params)

#EOF
