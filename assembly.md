---
---

# dDocent Reference Assembly Tutorial
Designed by Jon Puritz

**NOTE: You can download the RefTut file from the repository and run this tutorial from the command line**

# GOALS
1. To demultiplex samples with process_radtags and rename samples 
2. To use the methods of dDocent (via rainbow) to assemble reference contigs
3. To learn how optimize a de novo reference assembly

# Tutorial
*dDocent must be properly installed for this tutorial to work*


## Tutorial
*dDocent must be properly installed for this tutorial to work*

Start by downloading a small test dataset

```
curl -L -o data.zip https://www.dropbox.com/s/t09xjuudev4de72/data.zip?dl=0
```

Let's check that everything went well.

```bash
unzip data.zip
ll
```

You should see something like this:

```
Archive:  data.zip
  inflating: SimRAD.barcodes         
  inflating: SimRAD_R1.fastq.gz      
  inflating: SimRAD_R2.fastq.gz      
  inflating: simRRLs2.py             
total 7664
-rw-rw-r--. 1 j.puritz j.puritz 3127907 Mar  1 10:26 data.zip
-rwxr--r--. 1 j.puritz j.puritz     600 Mar  6  2015 SimRAD.barcodes
-rwxr--r--. 1 j.puritz j.puritz 2574784 Mar  6  2015 SimRAD_R1.fastq.gz
-rwxr--r--. 1 j.puritz j.puritz 2124644 Mar  6  2015 SimRAD_R2.fastq.gz
-rwxr--r--. 1 j.puritz j.puritz   12272 Mar  6  2015 simRRLs2.py
```

The data that we are going to use was simulated using the simRRLs2.py script that I modified from the one published by Deren Eaton.  You can find the original version [here](http://dereneaton.com/software/simrrls/).  Basically, the script simulated ddRAD 1000 loci shared across an ancestral population and two extant populations.  Each population had 180,000 individuals, and the two extant 
population split from the ancestral population 576,000 generations ago and split from each other 288,000 generation ago.  The two populations exchanged 4N*0.001 migrants per generation until about 2,000 generations ago.  4Nu equaled 0.00504 and mutations had a 10% chance of being an INDEL polymorphism.  Finally, reads for each locus were simulated on a per individual basis at a mean of 20X coverage (coming from a normaldistribution with a SD 8) and had an inherent sequencing error rate of 0.001. 

In short, we have two highly polymorphic populations with only slight levels of divergence from each other.  GST should be approximately
0.005. The reads are contained in the two fastq.gz files.

Let's go ahead and demultiplex the data.  This means we are going to separate individuals by barcode.
My favorite software for this task is process_radtags from the Stacks package (http://creskolab.uoregon.edu/stacks/) process_radtags takes fastq or fastq.gz files as input along with a file that lists barcodes.  Data can be separated according to inline
barcodes (barcodes in the actual sequence), Illumina Index, or any combination of the two.  Check out the manual at this website (http://creskolab.uoregon.edu/stacks/comp/process_radtags.php)

Let's start by making a list of barcodes.  The SimRAD.barcodes file actually has the sample name and barcode listed.  See for yourself.


You should see:

```
PopA_01 ATGGGG
PopA_02 GGGTAA
PopA_03 AGGAAA
PopA_04 TTTAAG
PopA_05 GGTGTG
PopA_06 TGATGT
PopA_07 GGTTGT
PopA_08 ATAAGT
PopA_09 AAGATA
PopA_10 TGTGAG
```

We need to turn this into a list of barcodes.  We can do this easily with the cut command.

```
cut -f2 SimRAD.barcodes > barcodes
```

Now we have a list of just barcodes.  The cut command let's you select a column of text with the -f (field command).  We used -f2 to get the second column.  

```
head barcodes
```

Now we can run process_radtags

```
process_radtags -1 SimRAD_R1.fastq.gz -2 SimRAD_R2.fastq.gz -b barcodes -e ecoRI --renz_2 mspI -r -i gzfastq
```

The option -e specifies the 5' restriction site and `--renze_2` specifes the 3' restriction site.  `-i` states the format of the input 
sequences.The `-r` option tells the program to fix cut sites and barcodes that have up to 1-2 mutations in them.  This can be changed 
with the `--barcode_dist flag`.  

Once the program is completed.  Your output directory should have several files that look like: 
`sample_AAGAGG.1.fq.gz, sample_AAGAGG.2.fq.gz, sample_AAGAGG.rem.1.fq.gz, and sample_AAGAGG.rem.2.fq.gz`

The *.rem.*.fq.gz files would normally have files that fail process_radtags (bad barcode, ambitious cut sites), but we have 
simulated data and none of those bad reads.  We can delete.
```
rm *rem*
```
The individual files are currently only names by barcode sequence.  We can rename them in an easier convention using a simple bash script.
Download the "Rename_for_dDocent.sh" script from my github repository

```
curl -L -O https://github.com/jpuritz/dDocent/raw/master/Rename_for_dDocent.sh
```

Take a look at this simple script

```
cat Rename_for_dDocent.sh
```

Bash scripts are a wonderful tool to automate simple tasks.  This script begins with an If statement to see if a file was provided as input.  If the file is not it exits and says why.  The file it requires is a two column list with the sample name in the first column and sample barcode in the second column.  The script reads all the names into an array and all the barcodes into a second array, and then gets the length of both arrays.  It then iterates with a for loop the task of renaming the samples.  

Now run the script to rename your samples and take a look at the output

```bash
bash Rename_for_dDocent.sh SimRAD.barcodes
ls *.fq.gz
```
There should now be 40 individually labeled .F.fq.gz and 40 .R.fq.gz.  Twenty from PopA and Twenty from PopB.
Now we are ready to rock!

Let's start by examining how the dDocent pipeline assembles RAD data.

First, we are going to create a set of unique reads with counts for each individual

```bash
ls *.F.fq.gz > namelist
sed -i'' -e 's/.F.fq.gz//g' namelist
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'

cat namelist | parallel --no-notice -j 8 "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
cat namelist | parallel --no-notice -j 8 "zcat {}.R.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
cat namelist | parallel --no-notice -j 8 "paste -d '-' {}.forward {}.reverse | mawk '$AWK3' | sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"
```

The first four lines simply set shell variables for various bits of AWK and perl code, to make parallelization with GNU-parallel easier. The first line after the variables, creates a set of forward reads for each individual by using mawk (a faster, c++ version of awk) to sort through the fastq file and strip away the quality scores.  The second line does the same for the PE reads.  Lastly, the final line concatentates the forward and PE reads together (with 10 Ns between them) and then find the unique reads within that individual and counts the occurences (coverage).

Now we can take advantage of some of the properties for RAD sequencing.  Sequences with very small levels of coverage within an individual are likely to be sequencing errors.  So for assembly we can eliminate reads with low copy numbers to remove non-informative data!

Let's sum up the number the within individual coverage level of unique reads in our data set

```bash
cat *.uniq.seqs > uniq.seqs
for i in {2..20};
do 
echo $i >> pfile
done
cat pfile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data
rm pfile
```

This is another example of a BASH for loop.  It uses mawk to query the first column and
select data above a certain copy number (from 2-20) and prints that to a file.

Take a look at the contents of uniqseq.data

```bash
more uniqseq.data
```

We can even plot this to the terminal using gnuplot

```bash
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [2:20] 
unset label
set title "Number of Unique Sequences with More than X Coverage (Counted within individuals)"
set xlabel "Coverage"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.data' with lines notitle
pause -1
EOF
```

```bash
                       Number of Unique Sequences with More than X Coverage (Counted within individuals)
  Number of Unique Sequences
    70000 ++----------+-----------+-----------+-----------+----------+-----------+-----------+-----------+----------++
          +           +           +           +           +          +           +           +           +           +
          |                                                                                                          |
    60000 ******                                                                                                    ++
          |     ******                                                                                               |
          |           ******                                                                                         |
          |                 ******                                                                                   |
    50000 ++                      *****                                                                             ++
          |                            *                                                                             |
          |                             *****                                                                        |
    40000 ++                                 *                                                                      ++
          |                                   ******                                                                 |
          |                                         *****                                                            |
          |                                              *                                                           |
    30000 ++                                              *****                                                     ++
          |                                                    *                                                     |
          |                                                     *****                                                |
    20000 ++                                                         ******                                         ++
          |                                                                ******                                    |
          |                                                                      ******                              |
          |                                                                            ******                        |
    10000 ++                                                                                 ************           ++
          |                                                                                              *************
          +           +           +           +           +          +           +           +           +           +
        0 ++----------+-----------+-----------+-----------+----------+-----------+-----------+-----------+----------++
          2           4           6           8           10         12          14          16          18          20
                                                           Coverage
```

Now we need to choose a cutoff value.
We want to choose a value that captures as much of the diversity of the data as possible 
while simultaneously eliminating sequences that are likely errors.
Let's try 4

```bash
parallel --no-notice -j 8 mawk -v x=4 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv
wc -l uniqCperindv
```

We've now reduced the data to assemble down to 7598 sequences!
But, we can go even further.
Let's now restrict data by the number of different individuals a sequence appears within.

```bash
for ((i = 2; i <= 10; i++));
do
echo $i >> ufile
done

cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data
rm ufile
```

Again, we can plot the data:

```bash
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Number of Unique Sequences present in more than X Individuals"
set xlabel "Number of Individuals"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.peri.data' with lines notitle
pause -1
EOF
```

```bash
                                 Number of Unique Sequences present in more than X Individuals
  Number of Unique Sequences
    6000 ++------------+------------+-------------+------------+-------------+------------+-------------+-----------++
         +             +            +             +            +             +            +             +            +
         |                                                                                                           |
    5500 *****                                                                                                      ++
         |    **                                                                                                     |
    5000 ++     ***                                                                                                 ++
         |         ***                                                                                               |
         |            *                                                                                              |
    4500 ++            *****                                                                                        ++
         |                  ****                                                                                     |
         |                      ***                                                                                  |
    4000 ++                        *                                                                                ++
         |                          ***********                                                                      |
    3500 ++                                    ***                                                                  ++
         |                                        **********                                                         |
         |                                                  ***                                                      |
    3000 ++                                                    ***********                                          ++
         |                                                                ***                                        |
         |                                                                   *************                           |
    2500 ++                                                                               **************            ++
         |                                                                                              *************|
    2000 ++                                                                                                         +*
         |                                                                                                           |
         +             +            +             +            +             +            +             +            +
    1500 ++------------+------------+-------------+------------+-------------+------------+-------------+-----------++
         2             3            4             5            6             7            8             9            10
                                                     Number of Individuals

```

Again, we need to choose a cutoff value.
We want to choose a value that captures as much of the diversity of the data as possible 
while simultaneously eliminating sequences that have little value on the population scale.
Let's try 4.

```bash
mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.4.seqs
wc -l uniq.k.4.c.4.seqs
```

Now we have reduced the data down to only 3840 sequences!

Let's quickly convert these sequences back into fasta format
We can do this with two quick lines of code:

```bash
cut -f2 uniq.k.4.c.4.seqs > totaluniqseq
mawk '{c= c + 1; print ">Contig_" c "\n" $1}' totaluniqseq > uniq.fasta
```
This simple script reads the totaluniqseq file line by line and add a sequence header of >Contig X

### At this point, dDocent also checks for reads that have a substantial amount of Illumina adapter in them. Our data is simulated and does not contain adapter, so we'll skip that step for the time being.

With this, we have created our reduced data set and are ready to start assembling reference contigs.

First, let's extract the forward reads.

```bash
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
```

This uses the sed command to replace the 10N separator into a tab character and then uses the cut
function to split the files into forward reads.

Previous versions of dDocent utilized the program rainbow to do full RAD assembly; however, as of dDocent 2.0, parts of rainbow have been replaced for better functionality.  
For example, first step of rainbow clusters reads together using a spaced hash to estimate similarity in the forward reads only.  
dDocent now improves this by using clustering by alignment via the program CD-hit to achieve more accurate clustering.  Custom AWK code then converts the output of CD-hit to match the input of the 2nd phase of rainbow.

```bash
cd-hit-est -i uniq.F.fasta -o xxx -c 0.8 -T 0 -M 0 -g 1
```

This code clusters all off the forward reads by 80% similarity.  This might seem low, but other functions of rainbow will break up clusters given the number and frequency of variants, so it's best to use a low value at this step.

```bash
mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids
paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/\t/g' > rcluster
```

This code then converts the output of CD-hit to match the output of the first phase of rainbow.

The output follows a simple text format of:

```
Read_ID	Cluster_ID	Forward_Read	Reverse_Read
```

Use the more, head, and/or tail function to examine the output file (rcluster)
You should see approximately 1000 as the last cluster.  It's important to note that the numbers are not totally sequential and that there may not be 1000 clusters.  Try the command below to get the exact number.

```bash
cut -f2 rcluster | uniq | wc -l 
```

The actual number of clusters is 1000 in this case because this is simulated data.  

The next step of rainbow is to split clusters formed in the first step into smaller clusters representing significant variants.
Think of it in this way.  The first clustering steps found RAD loci, and this step is splitting the loci into alleles. 
This *also* helps to break up over clustered sequences.

```bash
rainbow div -i rcluster -o rbdiv.out 
```

The output of the div process is similar to the previous output with the exception that the second column is now the new divided cluster_ID
(this value is numbered sequentially) and there was a column added to the end of the file that holds the original first cluster ID
The parameter -f can be set to control what is the minimum frequency of an allele necessary to divide it into its own cluster
Since this is from multiple individuals, we want to lower this from the default of 0.2.

```
rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
```

Though changing the parameter for this data set has no effect, it can make a big difference when using real data.

The third part of the rainbow process is to used the paired end reads to merge divided clusters.  This helps to double check the clustering and dividing of the previous steps
all of which were based on the forward read.  The logic is that if divided clusters represent alleles from the same homolgous locus, they should have fairly similar paired end reads
as well as forward.  Divided clusters that do not share similarity in the paired-end read represent cluster paralogs or repetitive regions.  After the divided clusters are merged,
all the forward and reverse reads are pooled and assembled for that cluster.

```bash
rainbow merge -o rbasm.out -a -i rbdiv.out
```

A parameter of interest to add here is the -r parameter, which is the minimum number of reads to assemble.  The default is 5 which works well if assembling reads from a single individual.
However, we are assembling a reduced data set, so there may only be one copy of a locus.  Therefore, it's more appropriate to use a cutoff of 2.

```bash
rainbow merge -o rbasm.out -a -i rbdiv.out -r 2
```
The rbasm output lists optimal and suboptimal contigs.  Previous versions of dDocent used rainbow's included perl scripts to retrieve optimal contigs.  However, as of version 2.0, dDocent uses customized AWK code to extract optimal contigs for RAD sequencing.  

```bash
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
```

Now, this looks a bit complicated, but it's performing a fairly simple algorithm.  First, the script looks at all the contigs assembled for a cluster.  If any of the contigs contain forward and PE reads, then that contig is output as optimal.  If no overlap contigs exists (the usual for most RAD data sets), then the contig with the most assembled reads PE (most common) is output with the forward read contig with a 10 N spacer.  If two contigs have equal number of reads, the longer contig is output. 

### At this point, dDocent (version 2.0 and higher) will check for substantial overlap between F and PE reads in the contigs.  Basically double checking rainbow's assembly.  We will skip this for our simulated data though.###

Though rainbow is fairly accurate with assembly of RAD data, even with high levels of INDEL polymorphism.  It's not perfect and the resulting contigs need to be aligned
and clustered by sequence similarity.  We can use the program cd-hit to do this.

```bash
cd-hit-est -i rainbow.fasta -o referenceRC.fasta -M 0 -T 0 -c 0.9
```

The `-M` and `-T` flags instruct the program on memory usage (-M) and number of threads (-T).  Setting the value to 0 uses all available.  The real parameter of significan is the -c parameter which
sets the percentage of sequence similarity to group contigs by.  The above code uses 90%.  Try using 95%, 85%, 80%, and 99%.
Since this is simulated data, we know the real number of contigs, 1000.  By choosing an cutoffs of 4 and 4, we are able to get the real number of contigs, no matter what the similarty cutoff.  

In this example, it's easy to know the correct number of reference contigs, but with real data this is less obvious.  As you just demonstrated, varying the uniq sequence copy cutoff and the final clustering similarity have the
the largest effect on the number of final contigs.  You could go back and retype all the steps from above to explore the data, but scripting makes this easier.
I've made a simple bash script called remake_reference.sh that will automate the process.  

```bash
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/remake_reference.sh
```

You can remake a reference by calling the script along with a new cutoff value and similarity.

```bash
bash remake_reference.sh 4 4 0.90 PE 2
```

This command will remake the reference with a cutoff of 20 copies of a unique sequence to use for assembly and a final clustering value of 90%.
It will output the number of reference sequences and create a new, indexed reference with the given parameters.
The output from the code above should be "1000"
Experiment with some different values on your own.   
What you choose for a final number of contigs will be something of a judgement call.  However, we could try to heuristically search the parameter space to find an optimal value.
Download the script to automate this process

```bash
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/ReferenceOpt.sh
```

Take a look at the script ReferenceOpt.sh.  
This script uses  different loops to assemble references from an interval of cutoff values and c values from 0.8-0.98.  It take as a while to run, so I have pasted the output below for you.

```bash
#bash ReferenceOpt.sh 4 8 4 8 PE 16
```

```bash
                                          Histogram of number of reference contigs
  Number of Occurrences
    200 ++--------------+--------------+---------------+--------------+---------------+--------------+--------------++
        +               +              +               +      'plot.kopt.data' using (bin($1,binwidth)):(1.0)*********
    180 ++                                                                                           *              +*
        |                                                                                            *               *
        |                                                                                            *               *
    160 ++                                                                                           *              +*
        |                                                                                            *               *
    140 ++                                                                                           *              +*
        |                                                                                            *               *
        |                                                                                            *               *
    120 ++                                                                                           *              +*
        |                                                                                            *               *
    100 ++                                                                                           *              +*
        |                                                                                            *               *
     80 ++                                                                                           *              +*
        |                                                                                            *               *
        |                                                                                            *               *
     60 ++                                                                                           *              +*
        |                                             ************************************************               *
     40 ++                                            *                                              *              +*
        |                                             *                                              *               *
        |                                             *                                              *               *
     20 ++                                            *                                              *              +*
        ***********************************************+              +               +              *               *
      0 **************************************************************************************************************
       988             990            992             994            996             998            1000            1002
                                                 Number of reference contigs

Average contig number = 999.452
The top three most common number of contigs
X	Contig number
164	1000
19	998
18	999
The top three most common number of contigs (with values rounded)
X	Contig number
250	1000.0
```

You can see that the most common number of contigs across all iteration is 1000, but also that the top three occuring and the average are all within 1% of the true value
Again, this is simulated data and with real data, the number of exact reference contigs is unknown and you will ultimately have to make a judgement call.

Let's examine the reference a bit.

```bash
bash remake_reference.sh 4 4 0.90 PE 2
head reference.fasta
```

```
>dDocent_Contig_1
NAATTCCTCCGACATTGTCGGCTTTAAATAGCTCATAACTTGAGCCCAGGTAAAGACTTTAGTATACTCGCACCTTCCGCTTATCCCCCGGCCGCNNNNNNNNNNATTCAACCGCGGGACCTGAACTAACATAGCGTTGTGTATACCATCCGAGGTAACCTTATAACTCTCTGCCATTCGGACAGGTAACACGGCATATCGTCCGN
>dDocent_Contig_2
NAATTCAGAATGGTCATACAGGGCGGTAGAATGGAATCCTGAATCGAATGGCGGTTGCATTGAGAACCTGGTACCAGATAGGATCTGGATTAAATNNNNNNNNNNGTCGGGTACTAATTATCTATTGGGTCCAAACCCTCCGCCCCGTTTACTGCCCACCCGGCATGCAGTCATGAGAATTCCAAGGAACTAAGATAAGAGACCGN
>dDocent_Contig_3
NAATTCGGGCTCCTTGGAGAGATTCTTTCAATTATGCCCCCTACGTGGGAAACAGGGTCGGAAGTGGTCGGCTGAGAATTACTCGAAAGCCGCTCNNNNNNNNNNCCACCAGCATGATAGGACTTCAAGCTTGCCGTTTGTTGGGAGGACCGGTCGCTACGGAGCTGACGCTATCTCCCGCATCGGACCTCGTGGACAAAAACCGN
>dDocent_Contig_4
NAATTCAAAAGTCGCCCATAGGTACGTGATGAATTAGGTCAAGCGGGGACGTCGCATAGATGCGTGACGTCTGGAGCATGATGTTGTTTCTAACCNNNNNNNNNNAATCACTCGGTCAACGTGGTCCGTGCTCTGCAACGAAAAAAACTTCGCATGTGAACGATGATGCCTATAGGTGCGACCGCCGTCAGAGGCCCGTTGACCGN
>dDocent_Contig_5
NAATTCATACGGATATGATACTTCGTCTGGCAGGGTGGCTAGCGAGTTTAAGGATTCTTGGATAAAGGTAGGTAAAATTCTCGAGATTCTGATCTNNNNNNNNNNTAGAGGTGCTGGCGGGGCCTAGACGTGTTTCTACGCTTACTGATCAAATTAGCTAGCTTAGGTTCCTATAGTCTACGCTGGATTGTCCTTAGATGCACCGN
```
You can now see that we have complete RAD fragments starting with our EcoRI restriction site (AATT), followed by R1, then a filler of 10Ns,
and then R2 ending with the mspI restriction site (CCG). The start and end of the sequence are buffered with a single N

We can use simple shell commands to query this data.
Find out how many lines in the file (this is double the number of sequences)

```bash
wc -l reference.fasta
```

Find out how many sequences there are directly by counting lines that only start with the header character ">"

```bash
mawk '/>/' reference.fasta | wc -l 
```

We can test that all sequences follow the expected format.

```bash
mawk '/^NAATT.*N*.*CCGN$/' reference.fasta | wc -l
grep '^NAATT.*N*.*CCGN$' reference.fasta | wc -l
```

No surprises here from our simulated data, butI highly recommend familiarizing yourself with grep, awk, and regular expressions to help evaluate de novo references.

#Bonus Section

Here, I am going to let you in on an experimental script I have been using to help optimize reference assemblies.

```bash
curl -L -O https://raw.githubusercontent.com/jpuritz/WinterSchool.2016/master/RefMapOpt.sh
```

This script assembles references across cutoff values and then maps 20 random samples and evaluates mappings to the reference, along with number of contigs and coverage.  
It takes a long time to run, but here's a sample command and output

```bash
#RefMapOpt.sh 4 8 4 8 0.9 64 PE
```

This would loop across cutoffs of 4-8 using a similarity of 90% for clustering, parellized across 64 processors, using PE assembly technique.

The output is stored in a file called `mapping.results`

```bash
curl -L -o mapping.results https://www.dropbox.com/s/x7p7j1xn1hjltzv/mapping.results?dl=0
cat mapping.results
```

```
Cov		Non0Cov	Contigs	MeanContigsMapped	K1	K2	SUM Mapped	SUM Properly	Mean Mapped	Mean Properly	MisMatched
37.3382	39.6684	1000	942.25				4	4	747510		747342			37375.5		37367.1			0
37.4003	39.7343	1000	942.25				4	5	748753		748546			37437.7		37427.3			0
37.4625	39.7919	1000	942.45				4	6	749999		749874			37499.9		37493.7			0
37.4967	39.8282	1000	942.45				4	7	750685		750541			37534.2		37527.1			0
37.486	39.8169	1000	942.45				4	8	750469		750205			37523.4		37510.2			0
37.3517	39.6785	1000	942.35				5	4	747780		747612			37389		37380.6			0
37.4147	39.7454	1000	942.35				5	5	749042		748835			37452.1		37441.8			0
37.4701	39.7999	1000	942.45				5	6	750151		750009			37507.6		37500.4			0
37.4852	39.8161	1000	942.45				5	7	750453		750210			37522.7		37510.5			0
37.4551	39.7824	999		941.55				5	8	749102		748837			37455.1		37441.8			0
37.3561	39.6833	1000	942.35				6	4	747870		747731			37393.5		37386.6			0
37.453	39.7776	1000	942.55				6	5	749809		749734			37490.4		37486.7			0
37.4923	39.8193	1000	942.55				6	6	750595		750376			37529.8		37518.8			0
37.4784	39.8089	1000	942.45				6	7	750318		750075			37515.9		37503.8			0
37.4437	39.766	999		941.65				6	8	748874		748616			37443.7		37430.8			0
37.4013	39.7312	1000	942.35				7	4	748774		748698			37438.7		37434.9			0
37.4592	39.7907	1000	942.4				7	5	749934		749835			37496.7		37491.8			0
37.4682	39.7981	1000	942.45				7	6	750114		749897			37505.7		37494.8			0
37.4239	39.7468	1000	942.55				7	7	749227		748993			37461.3		37449.7			0
37.417	39.736	998		940.75				7	8	747591		747320			37379.6		37366			0
37.4413	39.761	1000	942.65				8	4	749575		749499			37478.8		37474.9			0
37.4492	39.7843	1000	942.3				8	5	749733		749562			37486.7		37478.1			0
37.4441	39.7711	998		940.6				8	6	748133		747888			37406.7		37394.4			0
37.4274	39.7517	997		939.7				8	7	747052		746779			37352.6		37338.9			0
37.5014	39.8269	989		932.25				8	8	742528		742279			37126.4		37113.9			0
```

I have added extra tabs for readability.  The output contains the average coverage per contig, the average coverage per contig not counting zero coverage contigs, the number of contigs, the mean number of contigs mapped, the two cutoff values used, the sum of all mapped reads, the sum of all properly mapped reads, the mean number of mapped reads, the mean number of properly mapped reads, and the number of reads that are mapped to mismatching contigs.
Here, we are looking to values that maximize properly mapped reads, the mean number of contigs mapped, and the coverage.  In this example, it's easy.  Values 4,7 produce the highes number of properly mapped reads, coverage, and contigs.  
Real data will involve a judgement call.  Again, I haven't finished vetting this script, so use at your own risk.

Congrats!  You've finished the reference assembly tutorial for dDocent.


