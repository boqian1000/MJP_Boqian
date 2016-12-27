for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_omh.py -n $i") > jo2.$i.sh; done
for i in jo2*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_mh.py -n $i") > jmh2.$i.sh; done
for i in jmh2*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_gbs.py -n $i") > jgbs2.$i.sh; done
for i in jgbs2*.sh; do echo $i; qsub $i; done
