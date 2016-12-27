for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_omh.py -n $i") > jo3.$i.sh; done
for i in jo3*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_mh.py -n $i") > jmh3.$i.sh; done
for i in jmh3*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_gbs.py -n $i") > jgbs3.$i.sh; done
for i in jgbs3*.sh; do echo $i; qsub $i; done
