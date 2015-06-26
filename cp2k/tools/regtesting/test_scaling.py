#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett

import os
import re
import numpy as np
import sys
from os import path
from subprocess import call


#===============================================================================
def main():
    if(len(sys.argv) < 4):
        print("Usage: test_scaling.py [--dry-run] <threshold> <exe_dir> <inp_file1> ... <inp_fileN>")
        sys.exit(1)

    if(sys.argv[1] == "--dry-run"):
        dry_run = True
        threshold = float(sys.argv[2])
        exe_dir = sys.argv[3]
        inp_files = sys.argv[4:]
    else:
        dry_run = False
        threshold = float(sys.argv[1])
        exe_dir = sys.argv[2]
        inp_files = sys.argv[3:]

    configs = [(8, 1, "popt"),
               (8, 1, "psmp"),
               (4, 2, "psmp"),
               (2, 4, "psmp"),
               (1, 8, "psmp"),
               (1, 8, "ssmp")]

    error = False
    max_diff_rel = 0.0
    for inp_fn in inp_files:
        print("\nRunning tests for %s ..."%inp_fn)
        runtimes = []
        summary = []
        for c in configs:
            label = "%s-%dx%d"%(c[2], c[0], c[1])
            out_fn = inp_fn.replace(".inp", "-%s.out"%label)
            try:
                if(not dry_run):
                    run_cp2k(exe_dir, inp_fn, out_fn, *c)
                t = parse_output(out_fn)
                runtimes.append(t)
                summary.append("%s: runtime %.f sec"%(label, t))
            except Exception as e:
                summary.append("%s: Error: %s"%(label, str(sys.exc_info()[1])))
                error = True

        print("\nTimings for "+inp_fn)
        print("  "+"\n  ".join(summary))
        mean = np.mean(runtimes)
        diff = np.max(np.abs(runtimes-mean))
        diff_rel = diff/mean * 100.0
        max_diff_rel = max(max_diff_rel, diff_rel)
        print("  Mean runtime: %.1f sec"%mean)
        print("  Max deviation: %.1f sec (%.1f%%)"%(diff, diff_rel))

    if(error):
        print("\nSummary: Something went wrong")
        print("Status: FAILED")
    else:
        print("\nSummary: Runtime varied by %.1f%% (threshold: %.1f%%)"%(max_diff_rel, threshold))
        print("Status: OK" if(max_diff_rel<threshold) else "Status: FAILED")


#===============================================================================
def run_cp2k(exe_dir, inp_fn, out_fn, nranks, nthreads, version):

    if(path.exists(out_fn)):
        os.remove(out_fn)
    launcher = "env OMP_NUM_THREADS=%d mpiexec -np %d %scp2k.%s"%(nthreads, nranks, exe_dir, version)
    cmd = launcher + " -i %s -o %s"%(inp_fn, out_fn)
    print("Running "+cmd)
    rtncode = call(cmd.split())
    if(rtncode != 0):
        print("Command exited with returncode: %d"%rtncode)
        raise(Exception("Run crashed"))

#===============================================================================
def parse_output(out_fn):
    output = open(out_fn).read()
    timing_report = re.search(r"\n( -+\n - +-\n - +T I M I N G +-\n([^\n]*\n){4}.* -+)\n", output, re.DOTALL).group(1)
    print timing_report+"\n"
    runtime = float(re.findall("\n CP2K      (.*)\n", output)[-1].split()[-1])
    return(runtime)

main()
#EOF
