for i in `seq 1 100`; do echo $i; (cat sub_hmc.sh; echo "python simu_server_hmc.py -n $i") > h1.$i.sh; done
for i in h1*.sh; do echo $i; qsub $i; done
