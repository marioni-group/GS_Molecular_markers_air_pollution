#GMRMomi for 8 pollutants & smoking - methylation data pre-regressed for variables

# Parallel BayesR 



#________________________________________________________
# Parallelise Step 1 
cd /gmrm-omi/build_gcc_openmpi
micromamba activate BayesR 


#Run this first
export OMP_NUM_THREADS=1


#Parallelising function
function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do #counts no background jobs, compares to no. specified when function called
wait -n #waits for one to finish if no. is >= no. specified
done
}


#Run
for i in /input_1yr/*.phen; 
do
A=$(basename -- "$i" .phen)
echo $i
echo $A

/home/jrober38/gmrm-omi/build_gcc_openmpi/gmrm \
--bin-file /input_1yr/full_residsc_prot15314.bin \
--dim-file /input_1yr/full_dim.dim \
--phen-files $i \
--group-index-file /input_1yr/full_gri.gri \
--group-mixture-file /input_1yr/grm.grm \
--shuffle-markers 1 \
--seed 171014 \
--iterations 2000 \
--out-dir /output_1yr/${A}/ &
pwait 2
done
wait
echo "All done" 



#Step 2________________________________________________________________________
cd /gmrm-omi/build_gcc_openmpi
micromamba activate BayesR

#Run this first
export OMP_NUM_THREADS=1

function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do
wait -n
done
}


for i in /input_1yr/*.phen;
do 
A=$(basename -- "$i" .phen)
echo $i
echo $A

/home/jrober38/gmrm-omi/build_gcc_openmpi/gmrm \
--bin-file /input_1yr/full_residsc_prot15314.bin  \
--dim-file /input_1yr/full_dim.dim \
--phen-files $i \
--iterations 2000 \
--burn-in 750 \
--model linear \
--test \
--in-name-base $A \
--out-dir /output_1yr/${A} &

pwait 2
done
wait
echo "All done" 


#GMRMomi for 8 pollutants - protein data pre-regressed for variables

# Parallel BayesR 
####### 6 month data #########


#________________________________________________________
# Parallelise Step 1 
cd /gmrm-omi/build_gcc_openmpi
micromamba activate BayesR #this version re-installed 17/06/25


#Run this first
export OMP_NUM_THREADS=1


#Parallelising function
function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do #counts no background jobs, compares to no. specified when function called
wait -n #waits for one to finish if no. is >= no. specified
done
}


#Run
for i in /input_6month/*.phen; 
do
A=$(basename -- "$i" .phen)
echo $i
echo $A

/home/jrober38/gmrm-omi/build_gcc_openmpi/gmrm \
--bin-file /full_residsc_prot15314.bin \
--dim-file /full_dim.dim \
--phen-files $i \
--group-index-file /full_gri.gri \
--group-mixture-file /grm.grm \
--shuffle-markers 1 \
--seed 171014 \
--iterations 2000 \
--out-dir /output_6month/${A}/ &
pwait 2
done
wait
echo "All done" 



#Step 2________________________________________________________________________
cd /gmrm-omi/build_gcc_openmpi
micromamba activate BayesR

#Run this first
export OMP_NUM_THREADS=1

function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do
wait -n
done
}


for i in /input_6month/*.phen;
do 
A=$(basename -- "$i" .phen)
echo $i
echo $A

/home/jrober38/gmrm-omi/build_gcc_openmpi/gmrm \
--bin-file /input_1yr/full_residsc_prot15314.bin  \
--dim-file /input_1yr/full_dim.dim \
--phen-files $i \
--iterations 2000 \
--burn-in 750 \
--model linear \
--test \
--in-name-base $A \
--out-dir /output_6month/${A} &

pwait 2
done
wait
echo "All done" 


# Parallel BayesR 
####### 7 year data #########


#________________________________________________________
# Parallelise Step 1 
cd /gmrm-omi/build_gcc_openmpi
micromamba activate BayesR #this version re-installed 17/06/25


#Run this first
export OMP_NUM_THREADS=1


#Parallelising function
function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do #counts no background jobs, compares to no. specified when function called
wait -n #waits for one to finish if no. is >= no. specified
done
}


#Run
for i in /input_7yr/*.phen; 
do
A=$(basename -- "$i" .phen)
echo $i
echo $A

/home/jrober38/gmrm-omi/build_gcc_openmpi/gmrm \
--bin-file /full_residsc_prot15314.bin \
--dim-file /full_dim.dim \
--phen-files $i \
--group-index-file /input_1yr/full_gri.gri \
--group-mixture-file /grm.grm \
--shuffle-markers 1 \
--seed 171014 \
--iterations 2000 \
--out-dir /output_7yr/${A}/ &
pwait 2
done
wait
echo "All done" 



#Step 2________________________________________________________________________
cd /gmrm-omi/build_gcc_openmpi
micromamba activate BayesR

#Run this first
export OMP_NUM_THREADS=1

function pwait() {
while [ $(jobs -p | wc -l) -ge $1 ]; do
wait -n
done
}


for i in /input_7yr/*.phen;
do 
A=$(basename -- "$i" .phen)
echo $i
echo $A

/home/jrober38/gmrm-omi/build_gcc_openmpi/gmrm \
--bin-file /full_residsc_prot15314.bin  \
--dim-file /full_dim.dim \
--phen-files $i \
--iterations 2000 \
--burn-in 750 \
--model linear \
--test \
--in-name-base $A \
--out-dir /output_7yr/${A} &

pwait 2
done
wait
echo "All done" 



