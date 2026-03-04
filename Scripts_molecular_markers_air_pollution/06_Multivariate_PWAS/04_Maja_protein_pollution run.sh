#Maja protein pollution run

#MAJA Run


eval "$(micromamba shell hook --shell bash)"

#or micromamba activate

#4. Activate the environment
micromamba activate MAJA


screen 
source /maja-python/bin/activate

cd /home_dir

#Run cluster saved files: Scaled pollution
mpiexec -n 4 python -m mpi4py /MAJA/maja.py \
--n 15314 \
--p 133 \
--q 9 \
--iters 2000 \
--burnin 1000 \
--x /protein_covcorr_0401.zarr \
--y /poll15314_sc.txt \
--dir /output_0425 \
--diagnostics True \
--g 133


