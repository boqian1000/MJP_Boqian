for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_server_pmcmc.py -n $i") > p2.$i.sh; done
for i in p2*.sh; do echo $i; qsub $i; done
