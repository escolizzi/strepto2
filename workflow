#for time plots: use plot_competition_popsize.py

### for heatmaps ###

#reduce file size by keeping only last season (this takes a while)
for fi in strpt_competition__krep_0.*.txt; do lastline=$(tail -1 $fi); IFS=' ' read -r -a larray <<< "$lastline"; awk -v a=${larray[0]} '{if($1==a){print $0}}' $fi > ${fi%.txt}_short.txt; done

#move files from the *.7 into *.6 directory:
for fi in strpt_competition__krep_0.*short.txt; do cp $fi ../testgenome/${fi%_short.txt}_3_short.txt; done

#make heatmap:
python3 ~/origins_postdoc/strepto_project/strepto2/plot_competition_popsize_heatmap.py <plotname> <season duration (2500)> <nr of datafiles for point (3)> datafile1.1 datafile 1.2 ...

#get ab production rate from seasonal report
awk 'BEGIN{sporsum=0; absum=0; count=0; spo1=0; ab1=0 } {if($1>2500 && $1 % 2500 == 900){spo1=$2; ab1=$3} else if($1>2500 && $1 % 2500 == 0){sporsum=sporsum+$2; absum=absum+($3-ab1); count=count+1; ab1=0; spo1=0}} END{print sporsum, absum, count, absum/(sporsum*count*1600*317)}' strpt_competition__gen3.txt >>proddata2.txt
