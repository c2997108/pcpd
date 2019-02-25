#!/bin/bash

if [ "$2" = "" ]; then echo $0 <an input fastq file> <a reference fasta file>; exit 1; fi

input=$1
ref=$2

nunagi=""
minlen=250
maxlen=295
minq=30

function exit1() {
 echo "error"
 exit
}

function ddcapipejudge() {
 for i in ${PIPESTATUS[*]}; do
  if [ $i != 0 ]; then
   if [ $i != 141 ]; then
    exit1;
   fi;
  fi;
 done
}

export TMPDIR=/tmp
mkdir OUTPUT_0001_0001

echo "cutadapt"
echo "cutadapt" > OUTPUT_0001_0001/ddca_summary.txt
cutadapt -b file:"`dirname $0`"/adapter.fasta -o file_cut.fq $input -q $minq,$minq -m $minlen >> OUTPUT_0001_0001/ddca_summary.txt || exit1

echo "fastq_quality_trimmer"
echo -e "\nfastq_quality_trimmer" >> OUTPUT_0001_0001/ddca_summary.txt
fastq_quality_trimmer -Q 33 -t $minq -l $minlen -i file_cut.fq -o file_cut_trim.fq -v >> OUTPUT_0001_0001/ddca_summary.txt || exit1

echo "fastq_quality_filter"
echo -e "\nfastq_quality_filter" >> OUTPUT_0001_0001/ddca_summary.txt
fastq_quality_filter -q $minq -p 100 -i file_cut_trim.fq -o file_cut_trim_filt.fq -Q 33 -v >> OUTPUT_0001_0001/ddca_summary.txt || exit1

echo "fastx_trimmer"
echo -e "\nfastx_trimmer" >> OUTPUT_0001_0001/ddca_summary.txt
fastx_trimmer -Q 33 -i file_cut_trim_filt.fq -o file_cut_trim_filt_trim.fq -v -l $maxlen >> OUTPUT_0001_0001/ddca_summary.txt || exit1

echo fastq_to_fasta
fastq_to_fasta -Q33 -i file_cut_trim_filt_trim.fq > file_cut_trim_filt_trim.fasta || exit1

echo sort fasta
fasta_to_tab.pl file_cut_trim_filt_trim.fasta |awk -F'\t' '{print length($2)"\t"$1"\t"$2}'| \
                            sort -k1,1nr|awk -F'\t' '{print ">"$2; print $3}' > file_cut_trim_filt_trim_sort.fasta

echo cd-hit
cd-hit-est -i file_cut_trim_filt_trim_sort.fasta -o file_cut_trim_filt_trim_sort_cdhit.fasta -d 0 -T 8 -M 5000 -c 1 || exit1

echo makeblastdb
sed 's/\r//g' $ref | awk '$0~"^>"{gsub(/\t/," ",$0); print $0} $0!~"^>"{gsub(/[-,._]/,"",$0); gsub(/[^acgtACGT]/,"N",$0); print $0}' > file_rename.fasta
makeblastdb -in file_rename.fasta -dbtype nucl -hash_index -out ref || exit1

echo blastn
blastn -db ref -query file_cut_trim_filt_trim_sort_cdhit.fasta -out file_cut_trim_filt_trim_sort_cdhit.fasta.blast -outfmt 6 -num_threads 8  || exit1

echo convert
cat file_cut_trim.fq|awk 'NR%4==1{ORS="\t"; print $1} NR%4==2{ORS="\n"; print $1}'|sed 's/^@//' > file_cut_trim.fq.tab
cat file_cut_trim_filt_trim_sort_cdhit.fasta|awk 'NR%2==1{ORS="\t"; print $1} NR%2==0{ORS="\n"; print $1}'|sed 's/^>//' > file_cut_trim_filt_trim_sort_cdhit.fasta.tab

awk -F'\t' 'BEGIN{maxlen='$maxlen'; minlen='$minlen'}
            FILENAME==ARGV[1]{
             if(length($2)>=minlen){
              if(length($2)>maxlen){len=maxlen}else{len=length($2)};
              for(i=minlen;i<=len;i++){if(a[substr($2,1,i)]==""){a[substr($2,1,i)]=$1}else{a[substr($2,1,i)]=a[substr($2,1,i)]"###"$1} }
             }
            }
            FILENAME==ARGV[2]{
             $2=substr($2,1,maxlen); if(a[$2]==""){print $1"\t"$2"\t"}else if(a[$2]!~"###"){print $1"\t"$2"\t"a[$2]}else{print $1"\t"$2"\t"}
            }' file_cut_trim_filt_trim_sort_cdhit.fasta.tab file_cut_trim.fq.tab > file_cut_trim.fq.tab.cnt || exit1

cat file_cut_trim.fq.tab.cnt | awk -F'\t' '$3!=""{cnt[$3]++} END{for(i in cnt){print i"\t"cnt[i]}}' |sort -k2,2nr > file_cut_trim.fq.tab.cnt.tsv

cat file_cut_trim_filt_trim_sort_cdhit.fasta.blast|
                  awk -F'\t' '{
                    if(old =="" && oldid != $1){old = $0; oldid = $1; oldscr=$12}
                    else if(old !="" && oldid ==$1 && $12<oldscr){print old; old=""}
                    else if(old !="" && oldid ==$1 && $12>=oldscr){split(old,arr,"\t"); OFS="\t"; $2=arr[2]","$2; old=$0}
                    else if(old !="" && oldid !=$1){print old; old = $0; oldid = $1; oldscr=$12}
                   } END{if(old!=""){print old}}' > file_cut_trim_filt_trim_sort_cdhit.fasta.blast.ann


awk -F'\t' 'FILENAME==ARGV[1]{a[$1]=$2} FILENAME==ARGV[2]{print $0"\t"a[$1]}' file_cut_trim_filt_trim_sort_cdhit.fasta.blast.ann file_cut_trim.fq.tab.cnt.tsv > OUTPUT_0001_0001/file.tsv

if [ "$nunagi" = "" ];then nunagi=`cat OUTPUT_0001_0001/file.tsv|wc -l`; fi
if [ "$nunagi" -gt `cat OUTPUT_0001_0001/file.tsv|wc -l` ]; then nunagi=`cat OUTPUT_0001_0001/file.tsv|wc -l`; fi
echo -e "\n$nunagi" >> OUTPUT_0001_0001/ddca_summary.txt

echo "" >> OUTPUT_0001_0001/ddca_summary.txt
echo "numbers of species" >> OUTPUT_0001_0001/ddca_summary.txt
head -n $nunagi OUTPUT_0001_0001/file.tsv |cut -f 3|sort|uniq -c|sort -nr > temp
cat temp >> OUTPUT_0001_0001/ddca_summary.txt

echo "" >> OUTPUT_0001_0001/ddca_summary.txt
echo "percentage of each species" >> OUTPUT_0001_0001/ddca_summary.txt
awk '{nunagi='$nunagi'; $1=$1/nunagi*100; print $0}' temp >> OUTPUT_0001_0001/ddca_summary.txt


