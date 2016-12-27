for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_omh.py -n $i") > jo1.$i.sh; done
for i in jo1*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_mh.py -n $i") > jmh1.$i.sh; done
for i in jmh1*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_gbs.py -n $i") > jgbs1.$i.sh; done
for i in jgbs1*.sh; do echo $i; qsub $i; done
