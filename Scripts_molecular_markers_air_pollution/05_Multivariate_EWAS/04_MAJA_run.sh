#MAJA Run
#For full test runs/troubleshooting runs see scripts/pollution_analysis/maja/EWAS/full_run_2/04_MAJA_fullRun.sh

#in terminal
eval "$(micromamba shell hook --shell bash)"

#or micromamba activate
cd /home_dir

#4. Activate the environment - make sure in home folder first
micromamba activate MAJA


screen 
source /maja-python/bin/activate


#Run home saved files
# output folder = maja_ewas/output_020725, screen = maja3
mpiexec -n 4 python -m mpi4py /MAJA/maja_v2_samdef3.py \
--n 18512 \
--p 752722 \
--q 9 \
--iters 1500 \
--burnin 500 \
--x full_residsc_meth18512_T.zarr \
--y full_poll_ids_newphenorder_2905.txt \
--dir output_020725/ \
--diagnostics True \
--g 752722 \
--Vthres 1

##Ran with 40 skipped iterations