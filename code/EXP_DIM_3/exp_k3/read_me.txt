for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_server.py -n $i") > qsub3.$i.sh; done
for i in qsub3*.sh; do echo $i; qsub $i; done
