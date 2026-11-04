# SPECTRE - Parallel, Pseudo-Spectral Navier-Stokes Solver for Turbulent Channel Flow  

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
nano spectre_config.nml # Edit simulation parameters in this file
mpirun -np <num_procs> ./spectre
```

### Re-Execution

```bash
nano spectre_config.nml # Change parameters here. No need to recompile for new simulation parameters
make init
mpirun -np <num_procs> ./spectre
```

### NCI Gadi Execution

```bash
make
nano spectre_config.nml
qsub myNCIJob.sh
```
