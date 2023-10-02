for file in `ls MDM*_marker.xls`
do
	id=`basename $file _marker.xls`
	cat $file |sed '1d' |awk '$5 >0{print $0}'|sort -k2,2nr |head -50 |cut -f1 |tr '\n' '\t' |sed 's/\t$/\n/g' |sed "s/^/$id\tNA\t/g" >>MDM.gmt
done

for file in `ls MG*_marker.xls`
do
	id=`basename $file _marker.xls`
	cat $file |sed '1d' |awk '$5 >0{print $0}'|sort -k2,2nr |head -50 |cut -f1 |tr '\n' '\t' |sed 's/\t$/\n/g' |sed "s/^/$id\tNA\t/g" >>MG.gmt
done

cat MDM.gmt MG.gmt >MDM_vs_MG.gmt