# sort
Fastest known Sort4 and Sort3 algorithms for the x86-64 architecture

This repository contains the functions described in our article "Sorting Performance Improvements beyond Deepmind's AlphaDev". The article features how we built upon the work for AlphaDev, uncovering yet faster sorting routines for short vectors on the x86-64 architecture by adopting a different algorithmic approach, using a broader range of machine instructions, and leveraging machine learning based function discovery methods.

## Running the code

The source file is structured to minimize dependencies: Just compile the single cpp file with -O3 optimizations and run the resulting executable. It checks whether all sort functions do sort correctly, and then benchmarks the different functions by running each of them 2bn times.

## About the measurements

For each function, the elapsed time for sorting 500 vectors at 4 different offsets 1m times (in total 2bn sorts) is measured. For each of the 1m runs, the vectors are reinitialized. The reinitialization is included in the results and distorts the results slightly. We considered taking the measurements excluding the reinitializations, but found that the compiler (clang in our case) did not inline the sort functions to be benchmarked any more if chrono-functions are moved into the loop, and the overhead of the call incl. surrounding spill code distorting the results much more than the reinitialization. As the reinitialization overhead is the same for all sorting functions, it tends to underestimate the speedup obtained by the faster routines.

We consider the comparison with std::sort as unfair, as std::sort can handle vectors of any length. We included it for lack of a better reference, and welcome any suggestions for improvement.

We also welcome any further improvements to the routines as well as measurements for other microarchitectures than those mentioned in the article (Intel Sandy Bridge, Intel Coffee Lake, AMD Zen+, AMD Zen2).
