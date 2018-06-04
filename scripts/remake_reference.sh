export LC_ALL=en_US.UTF-8

if [[ -z "$5" ]]; then
echo "Usage is sh remake_reference.sh K1 K2 similarity% Assembly_Type Number_of_Processors"
exit 1
fi

if ! awk --version | fgrep -v GNU &>/dev/null; then
         awk=gawk
    else
         awk=awk
fi

if find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null| grep -q 'trim' ; then
	TRIMMOMATIC=$(find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null | head -1)
	else
    echo "The dependency trimmomatic is not installed or is not in your" '$PATH'"."
    NUMDEP=$((NUMDEP + 1))
	fi
	
if find ${PATH//:/ } -maxdepth 1 -name TruSeq2-PE.fa 2> /dev/null | grep -q 'Tru' ; then
	ADAPTERS=$(find ${PATH//:/ } -maxdepth 1 -name TruSeq2-PE.fa 2> /dev/null | head -1)
	else
    echo "The file listing adapters (included with trimmomatic) is not installed or is not in your" '$PATH'"."
    NUMDEP=$((NUMDEP + 1))
    fi

ATYPE=$4
NUMProc=$5
CUTOFF=$1
CUTOFF2=$2
simC=$3
NAMES=( `cat "namelist" `)

Reference(){

AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
AWK4='{for(i=0;i<$1;i++)print}'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
SED1='s/^[ \t]*//'
SED2='s/\s/\t/g'
CUTOFF=$1
CUTOFF2=$2
simC=$3
FRL=$(zcat ${NAMES[0]}.F.fq.gz | mawk '{ print length() | "sort -rn" }' | head -1)

if [ ${NAMES[@]:(-1)}.F.fq.gz -nt ${NAMES[@]:(-1)}.uniq.seqs ];then
	if [[ "$ATYPE" == "PE" || "$ATYPE" == "RPE" ]]; then
	#If PE assembly, creates a concatenated file of every unique for each individual in parallel
		cat namelist | parallel --no-notice -j $NUMProc "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
		cat namelist | parallel --no-notice -j $NUMProc "zcat {}.R.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
		if [ "$ATYPE" = "RPE" ]; then
			cat namelist | parallel --no-notice -j $NUMProc "paste {}.forward {}.reverse | sort -k1 -S 100M > {}.fr"
			cat namelist | parallel --no-notice -j $NUMProc "cut -f1 {}.fr | uniq -c > {}.f.uniq && cut -f2 {}.fr > {}.r"
			cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK4' {}.f.uniq > {}.f.uniq.e" 
			cat namelist | parallel --no-notice -j $NUMProc "paste -d '-' {}.f.uniq.e {}.r | mawk '$AWK3'| sed 's/-/NNNNNNNNNN/' | sed -e '$SED1' | sed -e '$SED2'> {}.uniq.seqs"
			rm *.f.uniq.e *.f.uniq *.r *.fr
		else
			cat namelist | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"
		fi
		rm *.forward
		rm *.reverse
	fi
	if [ "$ATYPE" == "SE" ]; then
	#if SE assembly, creates files of every unique read for each individual in parallel
		cat namelist | parallel --no-notice -j $NUMProc "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
	fi
	if [[ "$ATYPE" == "OL" || "$ATYPE" == "ROL" ]]; then
	#If OL assembly, dDocent assumes that the marjority of PE reads will overlap, so the software PEAR is used to merge paired reads into single reads
		for i in "${NAMES[@]}";
      do
      zcat $i.R.fq.gz | head -2 | tail -1 >> lengths.txt
      done	
    MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.txt| head -1)
		LENGTH=$(( $MaxLen / 3))
		for i in "${NAMES[@]}"
			do
			echo "Running PEAR on sample" $i
			pearRM -f $i.F.fq.gz -r $i.R.fq.gz -o $i -j $NUMProc -n $LENGTH &> $i.kopt.log
			done
		cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.assembled.fastq | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
	fi
	if [ "$ATYPE" == "HYB" ]; then
	#If HYB assembly, dDocent assumes some PE reads will overlap but that some will not, so the OL method performed and remaining reads are then put through PE method
		for i in "${NAMES[@]}";
      do
      zcat $i.R.fq.gz | head -2 | tail -1 >> lengths.txt
      done	
    MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.txt| head -1)
    LENGTH=$(( $MaxLen / 3))
		for i in "${NAMES[@]}"
			do
			pearRM -f $i.F.fq.gz -r $i.R.fq.gz -o $i -j $NUMProc -n $LENGTH &>kopt.log
			done
		cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.assembled.fastq | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
		
		cat namelist | parallel --no-notice -j $NUMProc "cat {}.unassembled.forward.fastq | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
		cat namelist | parallel --no-notice -j $NUMProc "cat {}.unassembled.reverse.fastq | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
		cat namelist | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.ua.seqs"
		rm *.forward
		rm *.reverse
	fi		
	
fi

#Create a data file with the number of unique sequences and the number of occurrences
if [[ "$ATYPE" == "RPE" || "$ATYPE" == "ROL" ]]; then
  parallel --no-notice -j $NUMProc mawk -v x=$CUTOFF \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | sed 's/NNNNNNNNNN/-/' >  total.uniqs
  cut -f 1 -d "-" total.uniqs > total.u.F
  cut -f 2 -d "-" total.uniqs > total.u.R
  paste total.u.F total.u.R | sort -k1 --parallel=$NUMProc -S 2G > total.fr
  special_uniq(){
    mawk -v x=$1 '$1 >= x' $2  |cut -f2 | sed -e 's/NNNNNNNNNN/	/g' | cut -f1 | uniq
  }
  export -f special_uniq
  parallel --no-notice --env special_uniq special_uniq $CUTOFF {} ::: *.uniq.seqs  | sort --parallel=$NUMProc -S 2G | uniq -c > total.f.uniq
  join -1 2 -2 1 -o 1.1,1.2,2.2 total.f.uniq total.fr | mawk '{print $1 "\t" $2 "NNNNNNNNNN" $3}' | mawk -v x=$CUTOFF2 '$1 >= x' > uniq.k.$CUTOFF.c.$CUTOFF2.seqs
  rm total.uniqs total.u.* total.fr total.f.uniq* 
  
else
	parallel --no-notice mawk -v x=$CUTOFF \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | mawk -v x=$CUTOFF '$1 >= x' > uniq.k.$CUTOFF.c.$CUTOFF2.seqs
fi
sort -k1 -r -n --parallel=$NUMProc -S 2G uniq.k.$CUTOFF.c.$CUTOFF.seqs |cut -f2 > totaluniqseq
mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' totaluniqseq > uniq.full.fasta
LENGTH=$(mawk '!/>/' uniq.full.fasta  | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
LENGTH=$(($LENGTH * 3 / 4))
$awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' uniq.full.fasta > uniq.fq
java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 uniq.fq uniq.fq1 ILLUMINACLIP:$ADAPTERS:2:30:10 MINLEN:$LENGTH &>/dev/null
mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' uniq.fq1 > uniq.fasta
mawk '!/>/' uniq.fasta > totaluniqseq
rm uniq.fq*

#If this is a PE assebmle
if [[ "$ATYPE" == "PE" || "$ATYPE" == "RPE" ]]; then
	#Reads are first clustered using only the Forward reads using CD-hit instead of rainbow
	if [ "$ATYPE" == "PE" ]; then
	  sed -e 's/NNNNNNNNNN/	/g' uniq.fasta | cut -f1 > uniq.F.fasta
	  CDHIT=$(python -c "print (max("$simC" - 0.1,0.8))")
	  cd-hit-est -i uniq.F.fasta -o xxx -c $CDHIT -T $NUMProc -M 0 -g 1 -d 100 &>cdhit.log
	  mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>dDocent_Contig_,...]//g' | sort -g -k1 -S 2G --parallel=$NUMProc > sort.contig.cluster.ids
	  paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
	  #CD-hit output is converted to rainbow format
	  sort -k2,2 -g contig.cluster.totaluniqseq -S 2G --parallel=$NUMProc | sed -e 's/NNNNNNNNNN/	/g' > rcluster
	  rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
	  rainbow merge -o rbasm.out -a -i rbdiv.out -r 2 -N10000 -R10000 -l 20 -f 0.75
	else
	  sed -e 's/NNNNNNNNNN/	/g' totaluniqseq | cut -f1 | sort --parallel=$NUMProc -S 2G| uniq | mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' > uniq.F.fasta
	  CDHIT=$(python -c "print (max("$simC" - 0.1,0.8))")
	  cd-hit-est -i uniq.F.fasta -o xxx -c $CDHIT -T $NUMProc -M 0 -g 1 -d 100 &>cdhit.log
  	mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>dDocent_Contig_,...]//g' | sort -g -k1 -S 2G --parallel=$NUMProc > sort.contig.cluster.ids
  	paste sort.contig.cluster.ids <(mawk '!/>/' uniq.F.fasta) > contig.cluster.Funiq
  	sed -e 's/NNNNNNNNNN/	/g' totaluniqseq | sort --parallel=$NUMProc -k1 -S 2G | mawk '{print $0 "\t" NR}'  > totaluniqseq.CN
  	join -t $'\t' -1 3 -2 1 contig.cluster.Funiq totaluniqseq.CN -o 2.3,1.2,2.1,2.2 > contig.cluster.totaluniqseq
	  #CD-hit output is converted to rainbow format
	  sort -k2,2 -g -S 2G --parallel=$NUMProc contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/	/g' > rcluster
	  rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
	  rainbow merge -o rbasm.out -a -i rbdiv.out -r 2 -N10000 -R10000 -l 20 -f 0.75
	fi
	#This AWK code replaces rainbow's contig selection perl script
  cat rbasm.out <(echo "E") |sed 's/[0-9]*:[0-9]*://g' | mawk ' {
		if (NR == 1) e=$2;
		else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
		else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
		else if ($1 ~/C/) clus=$2;
		else if ($1 ~/L/) len=$2;
		else if ($1 ~/S/) seq=$2;
		else if ($1 ~/N/) freq=$2;
		else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
		else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
		else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
		}' > rainbow.fasta

	seqtk seq -r rainbow.fasta > rainbow.RC.fasta
	mv rainbow.RC.fasta rainbow.fasta

	#The rainbow assembly is checked for overlap between newly assembled Forward and Reverse reads using the software PEAR
	sed -e 's/NNNNNNNNNN/	/g' rainbow.fasta | cut -f1 | $awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.F.fq
	sed -e 's/NNNNNNNNNN/	/g' rainbow.fasta | cut -f2 | $awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.R.fq

	seqtk seq -r ref.R.fq > ref.RC.fq
	mv ref.RC.fq ref.R.fq
	LENGTH=$(mawk '!/>/' rainbow.fasta | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
	LENGTH=$(( $LENGTH * 5 / 4))
	
	pearRM -f ref.F.fq -r ref.R.fq -o overlap -p 0.001 -j $NUMProc -n $LENGTH &>kopt.log

	rm ref.F.fq ref.R.fq

	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.assembled.fastq > overlap.fasta
	mawk '/>/' overlap.fasta > overlap.loci.names
	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.forward.fastq > other.F
	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.reverse.fastq > other.R
	paste other.F other.R | mawk '{if ($1 ~ />/) print $1; else print $0}' | sed 's/	/NNNNNNNNNN/g' > other.FR

	cat other.FR overlap.fasta > totalover.fasta

	rm *.F *.R
fi

if [[ "$ATYPE" == "HYB" ]];then
	parallel --no-notice mawk -v x=$CUTOFF \''$1 >= x'\' ::: *.uniq.ua.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | mawk -v x=$2 '$1 >= x' > uniq.k.$CUTOFF.c.$CUTOFF2.ua.seqs
	AS=$(cat uniq.k.$CUTOFF.c.$CUTOFF2.ua.seqs | wc -l)
	if [ "$AS" -gt 1 ]; then
		cut -f2 uniq.k.$CUTOFF.c.$CUTOFF2.ua.seqs > totaluniqseq.ua
		mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' totaluniqseq.ua > uniq.full.ua.fasta
		LENGTH=$(mawk '!/>/' uniq.full.ua.fasta  | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
		LENGTH=$(($LENGTH * 3 / 4))
		$awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' uniq.full.ua.fasta > uniq.ua.fq
		java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 uniq.ua.fq uniq.ua.fq1 ILLUMINACLIP:$ADAPTERS:2:30:10 MINLEN:$LENGTH &>/dev/null
		mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' uniq.ua.fq1 > uniq.ua.fasta
		mawk '!/>/' uniq.ua.fasta > totaluniqseq.ua
		rm uniq.ua.fq*
		#Reads are first clustered using only the Forward reads using CD-hit instead of rainbow
		sed -e 's/NNNNNNNNNN/	/g' uniq.ua.fasta | cut -f1 > uniq.F.ua.fasta
		CDHIT=$(python -c "print(max("$simC" - 0.1,0.8))")
		cd-hit-est -i uniq.F.ua.fasta -o xxx -c $CDHIT -T 0 -M 0 -g 1 -d 100 &>cdhit.log
		mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>dDocent_Contig_,...]//g' | sort -g -k1 -S 2G --parallel=$NUMProc > sort.contig.cluster.ids.ua
		paste sort.contig.cluster.ids.ua totaluniqseq.ua > contig.cluster.totaluniqseq.ua
		sort -k2,2 -g -S 2G --parallel=$NUMProc contig.cluster.totaluniqseq.ua | sed -e 's/NNNNNNNNNN/	/g' > rcluster.ua
		#CD-hit output is converted to rainbow format
		rainbow div -i rcluster.ua -o rbdiv.ua.out -f 0.5 -K 10
		if [ "$ATYPE" == "PE" ]; then
			rainbow merge -o rbasm.ua.out -a -i rbdiv.ua.out -r 2 -N10000 -R10000 -l 20 -f 0.75
		else
			rainbow merge -o rbasm.ua.out -a -i rbdiv.ua.out -r 2 -N10000 -R10000 -l 20 -f 0.75
		fi
		
		#This AWK code replaces rainbow's contig selection perl script
		cat rbasm.ua.out <(echo "E") |sed 's/[0-9]*:[0-9]*://g' | mawk ' {
			if (NR == 1) e=$2;
			else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_Contig_UA_" e "\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
			else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">dDocent_Contig_UA_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
			else if ($1 ~/C/) clus=$2;
			else if ($1 ~/L/) len=$2;
			else if ($1 ~/S/) seq=$2;
			else if ($1 ~/N/) freq=$2;
			else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
			else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
			else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
			}' > rainbow.ua.fasta
	
		seqtk seq -r rainbow.ua.fasta > rainbow.RC.fasta
		mv rainbow.RC.fasta rainbow.ua.fasta
	
		cat rainbow.ua.fasta uniq.fasta > totalover.fasta

	fi
fi

if [[ "$ATYPE" != "PE" && "$ATYPE" != "RPE" && "$ATYPE" != "HYB" ]]; then
	cp uniq.fasta totalover.fasta
fi

cd-hit-est -i totalover.fasta -o reference.fasta.original -M 0 -T $NUMProc -c $simC &>cdhit.log

sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.fasta.original > reference.fasta

SEQS=$(cat reference.fasta | wc -l)
SEQS=$(($SEQS / 2 ))
echo $SEQS
}
Reference $CUTOFF $CUTOFF2 $simC
