##########################Script to make DNAplotter graphics from .wig expression files###########

######Combining and comparing hfq and wt sequences##############
####Documents in Documents/Bioinformatics/DNAplotGraphics/
########################################################
##Prepair files by simplifying names or change script logic#####

###Trim and rename .wig files##
for i in *.wig; do tail -n +3 $i | awk '{print $1,$2}' > ${i}_count.txt; done

#######Join multiple files into same summery file according to .txt file with one nr per bp in genome (nrFile.txt)
#############If no such file make it with: 
seq -f %.0f 0 2272360 > nrFile.txt
#################################################################
##########reverse direction files##################
####join .wig#####
join -e0 -a1 -a1 -o0,2.2 nrFile.txt P5004_101_r.wig_count.txt | join -e 0 -a 1 -o0,1.2,2.2 - P5004_104_r.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,2.2 - P5004_107_r.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,2.2 - P5004_102_r.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,2.2 - P5004_105_r.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,1.6,2.2 -  P5004_108_r.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,1.6,1.7,2.2 - P5004_103_r.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2 - P5004_106_r.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2 - P5004_109_r.wig_count.txt > datareverse.txt
######calc mean#####
awk '{ print ($2 + $3 + $4)/3, ($5 + $6 + $7)/3 }' datareverse.txt > datareversemean.txt

############forward direction files##################
join -e0 -a1 -a1 -o0,2.2 nrFile.txt P5004_101_f.wig_count.txt | join -e 0 -a 1 -o0,1.2,2.2 - P5004_104_f.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,2.2 - P5004_107_f.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,2.2 - P5004_102_f.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,2.2 - P5004_105_f.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,1.6,2.2 -  P5004_108_f.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,1.6,1.7,2.2 - P5004_103_f.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2 - P5004_106_f.wig_count.txt | join -e 0 -a 1 -o0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2 - P5004_109_f.wig_count.txt > dataforward.txt
######calc mean#####
awk '{ print ($2 + $3 + $4)/3, ($5 + $6 + $7)/3 }' dataforward.txt > dataforwardmean.txt

#alt1##############differentiation calc through substraction######################## 
awk '{ print (($2 * -1) - ($1 * -1)) }' datareversemean.txt > datareversemeandifference.txt
awk '{ print ($2 - $1) }' dataforwardmean.txt > dataforwardmeandifference.txt
###Alt2
awk '{ if($1>$2) {print ($1 / $2} else {print ($2 / $1)} }' dataforwardmean.txt > dataforwardmeandifference.txt

#alt2##############differentiation calc through FC######################## 
awk '{ if (($1<(-1)) or ($2<(-1)) {print(0); next } if($1>$2) {print (($1 * -1) / ($2 * -1))} else {print (($2 * -1) / ($1 * -1))} }' datareversemean.txt > datareversemeandifference.txt
awk '{ if (($1>(-1)) or ($2>(-1))) {print(0); next } if($1>$2) {print (($1 * -1) / ($2 * -1))} else {print (($2 * -1) / ($1 * -1))} }' datareversemean.txt > datareversemeandifference.txt

######Combine read direction files into one######
paste dataforwardmeandifference.txt datareversemeandifference.txt > dataTotalMeanDifference.txt

#######add together strand values addition#######
awk '{ print ($1 + $2) }' dataTotalMeanDifference.txt > dataTotalMeanDifference_F.txt
#######add together strand values addition log#######
awk '{print (log($2 + 1))}' dataTotalMeanDifference.txt > dataTotalMeanDifferencelogHFQ.txt

#######Increase difference values by lowering min differenting expression, Needed becouse DNAplotter sucks#######
##Alt1
awk '{if ($1 < 0)  { print ((log($1 * -1 + 1)) * -1) } else {print (log($1 + 1))}}' dataTotalMeanDifference_F.txt > dataTotalMeanDifference_F_log.txt
###Alt2
awk '$1<0.2 {print($1 - 1); next} {print ($1 + 1)}' dataBothMeanDifference_FC_logged.txt > dataBothMeanDifference_FC_logged_And_Modded2.txt
