#03_MAJA_episcore_run


#set up


eval "$(micromamba shell hook --shell bash)"

#or micromamba activate
cd /home_dir

#4. Activate the environment - make sure in home folder first
micromamba activate MAJA


screen 
source /maja-python/bin/activate

# output folder = maja_episcore/output
mpiexec -n 4 python -m mpi4py /MAJA/maja_v2_samdef3.py \
--n 15516 \
--p 752722 \
--q 9 \
--iters 1500 \
--burnin 500 \
--x /full_residsc_meth15516T.zarr \
--y /poll_scaled_15516.txt \
--dir /maja_episcore/output/ \
--diagnostics True \
--g 752722 \
--Vthres 1

#Copy output to cluster
