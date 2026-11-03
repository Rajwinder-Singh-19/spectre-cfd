# SPECTRE - Parallel, Pseudo-Spectral, Turbulent, Incompressible Flow Navier-Stokes Solver  

```bash
cd spectre_bin # Change the directory before running all the commands below
```

## COMPILATION INSTRUCTIONS

### New Compilation

```bash
make 
```

### Re-Compilation

```bash
make clean
make 
```

## EXECUTION INSTRUCTIONS

### Fresh Execution After Compiling

```bash
nano config.nml # Edit simulation parameters in this file
mpirun -np <num_procs> ./main
```

### Re-Execution

```bash
nano config.nml # Change parameters here. No need to recompile for new simulation parameters
make init
mpirun -np <num_procs> ./main
```

### NCI Gadi Execution

```bash
make
nano config.nml
qsub myNCIJob.sh
```
