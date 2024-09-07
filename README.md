# MPI-Parallelized MD Code

In the `src/main.cpp` file, you will find a very simple molecular dynamics code that uses a Lennard-Jones potential to account for interactions between particles.
Your task is to parallelize this code.
Do some profiling work (using `MPI_Wtime`) to assist in identifying which regions of the code dominate the cost of the calculation, and to provide you with feedback on the effectiveness of your parallelization work.

Ensure that in your parallelized code, the amount of memory required for the computation does not substantially increase when the code is run on larger numbers of ranks.
In particular, no rank should ever hold the entire set of particle coordinates simultaneously (the same rule applies to the forces and velocities).
You'll need to think carefully about how you parallelize the double loop over particles (lines `181` and `182`).

The initial state of the system is generated in such a way that you should get the same results regardless of the number of ranks you run on (at least, within the span of short simulations, such as we are using here).
