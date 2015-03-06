#dDocent Reference Assembly Tutorial
Designed by Jon Puritz

**NOTE: You can download the RefTut file from the depository and run this tutorial from the command line**

#GOALS
1. To demultiplex samples with process_radtags and rename samples 
2. To use the methods of dDocent (via rainbow) to assemble reference contigs
3. To learn how optimize a de novo reference assembly

#Tutorial
*dDocent must be properly installed for this tutorial to work*

Start by downloading a small test dataset
```
curl -L -o data.zip https://www.dropbox.com/sh/nh2km7n2k8egmge/AABWKxCXww4BKMZIcQXV6Vxma?dl=1
```
Let's check that everything went well.
```bash
unzip data.zip && ll
```
You should see something like this:
```
     total 4612
    -rwxr--r--. 1 jpuritz users 600 Feb 23 03:36 SimRAD.barcodes
    -rwxr--r--. 1 jpuritz users 2574784 Feb 23 03:36 SimRAD_R1.fastq.gz
    -rwxr--r--. 1 jpuritz users 2124644 Feb 23 03:36 SimRAD_R2.fastq.gz
    -rwxr--r--. 1 jpuritz users 12272 Feb 23 03:36 simRRLs2.py
```
The data that we are going to use was simulated using the simRRLs2.py script that I modified from the one published by Deren Eaton.  You can find the original version here (http://dereneaton.com/software/simrrls/).  Basically, the script simulated ddRAD 1000 loci shared across an ancestral population and two extant populations.  Each population had 180,000 individuals, and the two extant 
population split from the ancestral population 576,000 generations ago and split from each other 288,000 generation ago.  The two populations exchanged 4N*0.001 migrants per generation until about 2,000 generations ago.  4Nu equaled 0.00504 and mutations had a 10% chance of being an INDEL polymorphism.  Finally, reads for each locus were simulated on a per individual basis at a mean of 20X coverage (coming from a normaldistribution with a SD 8) and had an inherent sequencing error rate of 0.001. 

In short, we have two highly polymorphic populations with only slight levels of divergence from each other.  GST should be approximately
0.005. The reads are contained in the two fastq.gz files.

Let's go ahead and demultiplex the data.  This means we are going to separate individuals by barcode.
My favorite software for this task is process_radtags from the Stacks package (http://creskolab.uoregon.edu/stacks/) process_radtags takes fastq or fastq.gz files as input along with a file that lists barcodes.  Data can be separated according to inline
barcodes (barcodes in the actual sequence), Illumina Index, or any combination of the two.  Check out the manual at this website (http://creskolab.uoregon.edu/stacks/comp/process_radtags.php)

Let's start by making a list of barcodes.  The SimRAD.barcodes file actually has the sample name and barcode listed.  See for yourself.
```
head SimRAD.barcodes
```
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
The option -e specifies the 5' restriction site and `--renze_2` species the 3' restriction site.  `-i` states the format of the input 
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
Bash scripts are a wonderful tool to automate simple tasks.  This script begins with an If statement to see if a file was provided as 
input.  If the file is not it exits and says why.  The file it requires is a two column list with the sample name in the first column 
and sample barcode in the second column.  The script reads all the names into an array and all the barcodes into a second array, and 
then gets the length of both arrays.  It then iterates with a for loop the task of renaming the samples.  

Now run the script to rename your samples and take a look at the output
```bash
sh Rename_for_dDocent.sh SimRAD.barcodes
ls *.fq.gz
```
There should now be 40 individually labeled .F.fq.gz and 40 .R.fq.gz.  Twenty from PopA and Twenty from PopB.
Now we are ready to rock!

Let's start by examining how the dDocent pipeline assembles RAD data using the program rainbow
You'll notice that no trimming or adapter cleaning is done at this point.  This is because rainbow requires sequences to be uniform in length.
Don't worry, we will still be able to remove a large amount from the data.
Let's start by pooling all the forward and reverse reads together and removing the quality scores
```bash
zcat *.F.fq.gz | mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' > forward
zcat *.R.fq.gz | mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' > reverse
```
This code uses mawk (a faster, c++ version of awk) to sort through the fastq file and remove only the quality scores but keep the header and
sequence

Next we can use the program seqtk to quickly reverse complement the pooled reverse reads to put everything in a 5' to 3' orientation
```
seqtk seq -r reverse > reverseRC
```
Now we can use a perl script to paste the forward and reverse reads together with forward reads with 10 Ns in between
```
mergefq.pl forward reverseRC concat.fasta
```
Let's look at the new merged sequences
```bash
head concat.fasta
```
```
>lane1_fakedata0_0 1:N:0:_1
AATTCGGCTTGCAACGCAAGTGACGATTCCCACGGACATAACTGATCTAAGTAACTTCCAAATCTGGGAATGGGATTTCATAATTAAGGACTATNNNNNNNNNNACGACGAGCAATCCACAGACCTAGGCCCATCGAAGCGTCTTATGATTGATAACATCAGAGGGGGATGGGAGGTCCTGCTGTCGCATGGGAGAATACACCG
>lane1_fakedata0_1 1:N:0:_1
AATTCGGCTTGCAACGCAAGTGACGATTCCCACGGACATAACTGATCTAAGTAACTTCCAAATCTGGGAATGGGATTTCATAATTAAGGACTATNNNNNNNNNNACGACGAGCAATCCACAGACCTAGGCCCATCGAAGCGTCTTATGATTGATAACATCAGAGGGGGATGGGAGGTCCTGCTGTCGCATGGGAGAATACACCG
>lane1_fakedata0_2 1:N:0:_1
AATTCGGCTTGCAACGCAAGTGACGATTCCCACGGACATAACTGATCTAAGTAACTTCCAAATCTGGGAATGGGATTTCATAATTAAGGACTATNNNNNNNNNNACGACGAGCAATCCACAGACCTAGGCCCATCGAAGCGTCTTATGATTGATAACATCAGAGGGGGATGGGAGGTCCTGCTGTCGCATGGGAGAATACACCG
>lane1_fakedata0_3 1:N:0:_1
AATTCGGCTTGCAACGCAAGTGACGATTCCCACGGACATAACTGATCTAAGTAACTTCCAAATCTGGGAATGGGATTTCATAATTAAGGACTATNNNNNNNNNNACGACGAGCAATCCACAGACCTAGGCCCATCGAAGCGTCTTATGATTGATAACATCAGAGGGGGATGGGAGGTCCTGCTGTCGCATGGGAGAATACACCG
>lane1_fakedata0_4 1:N:0:_1
AATTCGGCTTGCAACGCAAGTGACGATTCCCACGGACATAACTGATCTAAGTAACTTCCAAATCTGGGAATGGGATTTCATAATTAAGGACTATNNNNNNNNNNACGACGAGCAATCCACAGACCTAGGCCCATCGAAGCGTCTTATGATTGATAACATCAGAGGGGGATGGGAGGTCCTGCTGTCGCATGGGAGAATACACCG
```
You can now see that we have complete RAD fragments starting with our EcoRI restriction site (AATT), followed by R1, then a filler of 10Ns,
and then R2 ending with the mspI restriction site (CCG).

We can use simple shell commands to query this data.

Find out how many lines in the file (this is double the number of sequences)
```bash
wc -l concat.fasta
```
Find out how many sequences there are directly by counting lines that only start with the header character ">"
```bash
mawk '/>/' concat.fasta | wc -l 
```
We can test that all sequences follow the expected format.
```bash
mawk '/^AATT.*N*.*CCG$/' concat.fasta | wc -l
grep '^AATT.*N*.*CCG$' concat.fasta | wc -l
```
I highly recommend familiarizing yourself with grep, awk, and regular expressions!
Now, that we've checked the data we can continue with our assembly. 
At this point, we don't really need to keep track of the sequence names, so let's lose them.
    mawk '!/>/' concat.fasta > concat.seq
This mawk action basically gets all the lines from the fasta file that do NOT match the header character

Now let's find all the unique sequences in the data and how many times they occur
```bash
perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' concat.seq > uniq.seqs
```
Now we can take advantage of some of the properties for RAD sequencing.
Alleles that are shared among individuals are likely have a higher copy number in the data set.
Sequences with small copy number are likely to be sequencing errors.
So for assembly we can eliminate reads with low copy numbers to remove non-informative data!
Let's count up the number of unique reads per level of copy number
```bash
for ((i = 1; i <= 50; i++));
do
J=$(mawk -v x=$i '$1 >= x' uniq.seqs | wc -l)
echo -e "$i""\t""$J" >> uniqseq.data
done
```
This is another example of a BASH for loop.  It uses mawk to query the firt column and
select data above a certain copy number (from 1-50) and prints that to a file.

Take a look at the contents of uniqseq.data
```bash
more uniqseq.data
```
We can even plot this to the terminal using gnuplot
```bash
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Number of Unique Sequences with More than X Occurrences"
set xlabel "Number of Occurrences"
set ylabel "Number of Unique Sequences
plot 'uniqseq.data' with dots notitle
pause -1
EOF
```
```bash
                                     Number of Unique Sequences with More than X Occurrences
  Number of Unique Sequences
    140000 ++---------+---------+----------+---------+----------+----------+---------+----------+---------+---------++
           + .        +         +          +         +          +          +         +          +         +          +
           |                                                                                                         |
    120000 ++                                                                                                       ++
           |                                                                                                         |
           |                                                                                                         |
           |                                                                                                         |
    100000 ++                                                                                                       ++
           |                                                                                                         |
           |                                                                                                         |
     80000 ++                                                                                                       ++
           |                                                                                                         |
           |                                                                                                         |
           |                                                                                                         |
     60000 ++                                                                                                       ++
           |                                                                                                         |
           |                                                                                                         |
     40000 ++                                                                                                       ++
           |                                                                                                         |
           |                                                                                                         |
           |                                                                                                         |
     20000 ++                                                                                                       ++
           |   .                                                                                                     |
           +     . .  . . . . . . . .  . . . . . . . .  . . . . . . . . .  . . . . . . . .  . . . . . . . .  . . . . .
         0 ++---------+---------+----------+---------+----------+----------+---------+----------+---------+---------++
           0          5         10         15        20         25         30        35         40        45         50
                                                      Number of Occurrences
```
That's not terribly informative, other than sequencing error makes up a lot of the variation
in this data set.  We certainly don't want to use the reads that only appear once in the data set
for assembly.  Let's remove that line from the uniqseq.data file
```bash
sed -i -e "1d" uniqseq.data 
```
Let's look at that graph again
```bash
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Number of Unique Sequences with More than X Occurrences"
set xlabel "Number of Occurrences"
set ylabel "Number of Unique Sequences
plot 'uniqseq.data' with dots notitle
pause -1
EOF
```
```bash
                                    Number of Unique Sequences with More than X Occurrences
  Number of Unique Sequences
    13000 ++---------+---------+----------+----------+----------+---------+----------+----------+---------+---------++
          +          +         +          +          +          +         +          +          +         +          +
    12000 ++  .                                                                                                     ++
          |                                                                                                          |
          |                                                                                                          |
    11000 ++                                                                                                        ++
          |                                                                                                          |
    10000 ++                                                                                                        ++
          |                                                                                                          |
          |                                                                                                          |
     9000 ++                                                                                                        ++
          |                                                                                                          |
     8000 ++    .                                                                                                   ++
          |        . .                                                                                               |
     7000 ++           . .                                                                                          ++
          |                . .                                                                                       |
          |                    .  .                                                                                  |
     6000 ++                        . . .                                                                           ++
          |                               . . .                                                                      |
     5000 ++                                     . . .                                                              ++
          |                                            . . . .  .                                                    |
          |                                                       . . . . . .                                        |
     4000 ++                                                                  .  . . . . . .                        ++
          +          +         +          +          +          +         +          +       .  . . . . . . .  .     +
     3000 ++---------+---------+----------+----------+----------+---------+----------+----------+---------+------.-.+.
          0          5         10         15         20         25        30         35         40        45         50
                                                     Number of Occurrences
```
Now that looks like a much more informative graph.  Now we need to choose a cutoff value.
We want to choose a value that captures as much of the diversity of the data as possible 
while simultaneously eliminating sequences that have little value on the population scale.
Let's try 5
```bash
mawk -v x=5 '$1 >= x' uniq.seqs | cut -f 2 > totaluniqseq
wc -l totaluniqseq
```
We've now reduced the data to assemble down to 7481 sequences!

Let's quickly convert these sequences back into fasta format
We can use a BASH while loop to do this.
```bash
cat totaluniqseq | while read line
do
echo ">Contig"$i >>uniq.fasta
echo $line >> uniq.fasta
i=$(($i + 1))
done
```
This simple script reads the totaluniqseq file line by line and add a sequence header of >Contig X

Now let's split them back into forward and reverse reads
```bash
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f2 > uniq.R.fasta
```
This uses the sed command to replace the 10N separator into a tab character and then uses the cut
function to split the files into forward and reverse reads.

Now, let's revers complement the reverse reads to put them back into the original read orientation.
```bash
seqtk seq -r uniq.R.fasta > uniq.RC.fasta
rm uniq.R.fasta
```
With this, we have created our reduced data set and are ready to feed it into the program rainbow for assembly.
The first step of rainbow clusters reads together using a spaced hash to estimate similarity in the forward reads only.
```bash
rainbow cluster -1 uniq.F.fasta -2 uniq.RC.fasta > rcluster
```
We the output follows a simple text format of:
```
Read_ID	Cluster_ID	Forward_Read	Reverse_Read
```
Use the more, head, and/or tail function to examine the output file (rcluster)
You should see approximately 1252 as the last cluster.  It's important to note that the numbers are not totally sequential and that there
are not 1252 clusters.  Try the command below to get the exact number.
```bash
cut -f2 rcluster | uniq | wc -l 
```
The actual number of clusters is 1035.  

We can change the clustering threshold with the -m parameter, this is the maximum number of mismatches (default 4) in the spaced hash.
What happens if you change it to 6?
```bash
rainbow cluster -m 6 -1 uniq.F.fasta -2 uniq.RC.fasta > rcluster
```
Or 2?
```bash
rainbow cluster -m 2 -1 uniq.F.fasta -2 uniq.RC.fasta > rcluster
```
Try some other values on your own
You can see that the clustering is fairly insensitive to this parameter.  This is due to the level of polymorphism in the data set and the
accuracy and sensitivity of spaced-seed clustering.

The next step of rainbow is to split clusters formed in the first step into smaller clusters representing significant variants.
Think of it in this way.  The first clustering steps found RAD loci, and this step is splitting the loci into alleles. 
This *also* helps to break up over clustered sequences.
```bash
rainbow div -i rcluster -o rbdiv.out
```
The output of the div process is similar to the previous output with the exception that the second column is now the new divided cluster_ID
(this value is numbered sequentially) and there was a column added to the end of the file that holds the original first cluster ID
The parameter -f can be set to control what is the minimum frequency of an allele necessary to divide it into its own cluster
Since this is pooled data, we want to lower this from the default of 0.2.
```
rainbow div -i rcluster -o rbdiv.out -f 0.01
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
The rbasm output lists optimal and suboptimal contigs.  To extract only the optimal contigs, we can use an included perl script.
```bash
select_best_rbcontig_plus_read1.pl rbasm.out rbdiv.out >rainbow.fasta
```
This script has selected the longest contig for each cluster and checked if the forward and reverse reads overlap in that contig.
If they did not, per usual for RAD data, then it pastes the forward sequence with a 10 N spacer.

Though rainbow is fairly accurate with assembly of RAD data, even with high levels of INDEL polymorphism.  It's not perfect and the resulting contigs need to be aligned
and clustered by sequence similarity.  We can use the program cd-hit to do this.
```bash
cd-hit-est -i rainbow.fasta -o referenceRC.fasta -M 0 -T 0 -c 0.9
```
The `-M` and `-T` flags instruct the program on memory usage (-M) and number of threads (-T).  Setting the value to 0 uses all available.  The real parameter of significan is the -c parameter which
sets the percentage of sequence similarity to group contigs by.  The above code uses 90%.  Try using 95%, 85%, 80%, and 99%.
As you can see, the similarity parameter can have an effect on the number of final contigs, ranging from 1004 to 1035 in this example.
Since this is simulated data, we know the real number of contigs, 1000.  By choosing an original cutoff of 5 and manipulating other parameters (mainly -c of CD-hit), we were able get as close 0.04% away from the actual total.
On the other hand, we could have also been 3.5% away as well.  Not too shabby for a first pass.

In this example, it's easy to know the correct number of reference contigs, but with real data this is less obvious.  As you just demonstrated, varying the uniq sequence copy cutoff and the final clustering similarity have the
the largest effect on the number of final contigs.  You could go back and retype all the steps from above to explore the data, but scripting makes this easier.
I've made a simple bash script called remake_reference.sh that will automate the process.  
```bash
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/remake_reference.sh
```
You can remake a reference by calling the script along with a new cutoff value and similarity.
```bash
sh remake_reference.sh 13 0.85 
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
This script uses two different loops to assemble references from K cutoff 2-30 and c values from 0.8-0.98.  It take as a while to run, so I have pasted the output below for you.
I have pasted the output below though.
```bash
ReferenceOpt.sh 2 30
```
```bash
                                            Histogram of number of reference contigs
    Number of Occurrences
      50 ++------------+-------*******------------+-------------+-------------+-------------+------------+------------++
         +             +       *     *            +             'plot.kopt.data' using (bin($1,binwidth)):(1.0) ****** +
      45 ++                    *     *                                                                                ++
         |                     *     *                                                                                 |
         |                     *     *                                                                                 |
      40 ++                    *     *                                                                                ++
         |                     *     *                                                                                 |
      35 ++                    *     *******                                                                          ++
         |                     *     *     *************                                                               |
         |                     *     *     *     *     *                                                               |
      30 ++                    *     *     *     *     *                                                              ++
         |                     *     *     *     *     *                                                               |
      25 ++                    *     *     *     *     *                                                              ++
         |                     *     *     *     *     ******                                                          |
      20 ++               ******     *     *     *     *    *                                                         ++
         |          *******    *     *     *     *     *    *                                                          |
         |          *     *    *     *     *     *     *    *     *******                                              |
      15 ++         *     *    *     *     *     *     *    *******     *                                             ++
         |          *     *    *     *     *     *     *    *     *     *******                                        |
      10 ++         *     *    *     *     *     *     *    *     *     *     *******                                 ++
         |    *******     *    *     *     *     *     *    *     *     *     *     *                                  |
         ******     *     *    *     *     *     *     *    *     *     *     *     *                                  |
       5 *+   *     *     *    *     *     *     *     *    *     *     *     *     *                                 ++
         *    *     *  +  *    *     *     *     *+    *    *   + *     *     *     ************         +             +
       0 ********************************************************************************************************-----++
        980           990           1000         1010          1020          1030          1040         1050          1060
                                                   Number of reference contigs

Average contig number = 1007.49
The top three most common number of contigs
X	Contig number
15	1000
13	1004
12	999
The top three most common number of contigs (with values rounded)
X	Contig number
289 1000.0
1 	1100.0
```
You can see that the most common number of contigs across all iteration is 1000, but also that the top three occuring and the average are all within 1% of the true value
Again, this is simulated data and with real data, the number of exact reference contigs is unknown and you will ultimately have to make a judgement call.
Congrats!  You've finished the reference assembly tutorial.


