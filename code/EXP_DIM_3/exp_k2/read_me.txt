for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_server.py -n $i") > qsub2.$i.sh; done
for i in qsub2*.sh; do echo $i; qsub $i; done
