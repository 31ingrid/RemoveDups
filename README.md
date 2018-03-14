This R function is designed to help remove duplicated loci from a fasta file of selected loci to be used in analysis of RADseq data. The benefit of this method is that loci that were unintentionally split into two or duplicated can be removed.

The general protocol is designed to remove run through the STACKS pipeline (Catchen et al. 2013) with relaxed conditions (-m 5, default values). Then run populations with small values of -m, -r, and -p (e.g. 5, 0.5, and 0.5) because the goal is not to identify alleles, just whether loci exist or not. 

The cstacks module output provides a batch_x.catalog.tags.tsv. This file is essentially a list of loci identified for this dataset, plus other details. More details about the specific fields in that file are provided in the STACKS manual.

Once you have run through STACKs one time:

1. Use the Unix command ‘grep’ to create a new file based on the consensus sequences that will return one line for each locus:

grep consensus batch_x.catalog.tags.tsv>DATA_output

2. Manipulate that file (DATA_output) into fasta format and save it as a text file. There will be 1 line per locus.

(Here are a few details for that process. First, save it as a text file and open it with some text editor such as textmate. Then open it in excel and delete all files except the identification file (column #3) and the DNA sequence. Further manipulate it in excel to convert it to a fasta format, then convert to a .csv file.

3. Remove any loci that are not in haplotypes.tsv. 
When stacks generates the haplotypes.tsv file, it also removes any loci that were not present in the minimum number of populations specified by the (-p) parameter, and any other filtering you have asked for. By removing loci that are not in haplotypes.tsv, you are removing loci from the group you want to use as a baseline that do not meet your filtering requirements.

R code to generate a fasta file with only loci present in the haplotypes.tsv file:

Convert the haplotypes.tsv file to a .csv format in excel.

haps=read.csv("/path/to/data/batch_x.haplotypes.csv",header=TRUE)
fast1=read.csv("/path/to/data/DATA_output.csv",header=FALSE)
fast2=as.data.frame(fast1[haps$Catalog.ID,])
write.csv(fast2,file="/path/to/data/DATA_output_haps.txt")
Open the file DATA_ouput_haps.txt in a text editor like TextMate and change the suffix to .fa (DATA_output_haps.fa).

4. Use bowtie to test whether any of the loci are duplicates:

Create a bowtie database with bowtie-1.2 legacy:

./bowtie-build /path/to/data/ DATA_output_haps.fa /path/to/data/DATA_database1

Now align the database with the fasta file to look for duplicate sequences.

./bowtie -f -v 3 -p 1 -k 2 /path/to/data/DATA_database1 /path/to/data/ DATA_output_haps.fa DATA_alignments.txt

NOTE: v tells # of mismatches regardless of quality, -f fasta, -q FASTQ (.fq I think), -r is for raw data (.fq.gz although it was already process_radtags) -p number of processors, -k the number of valid alignments to report.

5. The file DATA_alignments.txt contains bowtie alignments; which loci match other loci. Most of the alignments are between loci and themselves, but there are matches. These duplicates need to be removed.  The function RemoveDups finds duplicates, as well as triplicates, and clusters of loci of any size that exist in the data. Then it removes all but one of the loci and saves it in a new .csv file that is the name of your initial fasta file with “NoDups” appended to the end. This file should not have a header (first row). Delete it if it is there.

The code is in R, and it is a function called RemoveDups. It takes the name of the bowtie output file (that has been converted to .csv format) and a fasta file with duplicates (that has been converted to .csv) and can be run as follows:

RemoveDups("DATA_alignments.csv","DATA_output_haps.csv")

These two example files (DATA_alignments.csv and DATA_output_haps.csv) have been included to use for testing. They are shortened versions of actual data files. Generally RADseq data will produce much larger datasets so the code may need to run for ~10 minutes.

Do not bother removing barcodes if you have them because they help aligning sequences later on.

The output will be called DATA_output_haps_noDups.csv. Convert it to DATA_output_haps_noDups.fa with a text editor.  Then this can be used to create a new bowtie database that can be used to set up a second round of stacks that has no repeated loci.

6. Make a new bowtie database
./bowtie-build /Volumes/LaCie/pollock_2nd/stacks/AK/DATA_output_haps_noDups.fa /path/to/data/DATA_database_noDups

7. The second round of the STACKS pipeline can start as follows, with 1 command like this per barcoded sample.

pstacks -p 8 -t bowtie -f ./samples/SAMPLENAME -o ./stacks -i 1 -m 5 –bound_high 0.05

References:
Catchen, J., Hohenlohe, P.A., Bassham, S., Amores, A. and Cresko, W.A., 2013. Stacks: an analysis tool set for population genomics. Molecular ecology, 22(11), pp.3124-3140.
