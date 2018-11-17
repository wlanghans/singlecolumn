a=("THL" "QT" "WTHL" "WQT" "THVD" "QTD" "THGRAD")
for i in ${a[@]}
do
ncl writemmcoefficients.ncl pathin='"/global/scratch/langhans/metamodel_mm3/"' casename='"CPBL2"' var=\"$i\"
done
a=("THL" "QT" "WTHL" "WQT" "THVD" "QTD" "QC")
for i in ${a[@]}
do
ncl writemmcoefficients.ncl pathin='"/global/scratch/langhans/metamodel_mm3/"' casename='"BOMEX"' var=\"$i\"
done
