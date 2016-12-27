for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_omh.py -n $i") > o3.$i.sh; done
for i in o3*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_mh.py -n $i") > mh3.$i.sh; done
for i in mh3*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_gbs.py -n $i") > gbs3.$i.sh; done
for i in gbs3*.sh; do echo $i; qsub $i; done
