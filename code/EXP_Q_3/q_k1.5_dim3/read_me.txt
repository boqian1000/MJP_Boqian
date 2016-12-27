for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_omh.py -n $i") > o1.$i.sh; done
for i in o1*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_mh.py -n $i") > mh1.$i.sh; done
for i in mh1*.sh; do echo $i; qsub $i; done


for i in `seq 1 100`; do echo $i; (cat sub_.sh; echo "python simu_gbs.py -n $i") > gbs1.$i.sh; done
for i in gbs1*.sh; do echo $i; qsub $i; done
