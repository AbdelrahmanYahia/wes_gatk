RGID=$(zcat /home/marc/Desktop/data/samples/complex_wes/father_ACC1-LIBR_L2_2.fq.gz | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
echo -e "$RGID"