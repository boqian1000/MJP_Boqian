for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_server_pmcmc.py -n $i") > p3.$i.sh; done
for i in p3*.sh; do echo $i; qsub $i; done
