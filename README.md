# COMP383
Computational Biology Pipeline Project

**Purpose**:
The provided python script named "Pipeline.py" should be capable of doing the following:
1. Retrieving transcriptomes from one patient donor 2 and 6 dpi given the SRA and converted to paired-end fastq files (This step is already achieved through the download_fastq method, which makes a url given the srr_number list on line 51).
2. Output the nuumber of reads in each transcriptome before and after bowtie2 mapping into the log file called PipelineProject.log found in the directory to be named "PipelineProject_<Aaron>_<Juco>", which is preemptively created for you in the script already in lines 37-48
3. Assemble the transcriptomes together to produce 1 assembly via SPAdes and a k-mer size of 99
4. Output the number of contigs with a length greater than 100 into the PipelineProject.log file as well as the length of the assembly
5. Retrieve the longest contig from the spades assembly using blastn and write the top 10 hits': subject accession (sacc), percent identity (pident), alignment length (length), start of alignment in query (qstart), end of alignment in query (qend), start of alignment in subject (sstart), end of alignment in subject (send), bitscore, e-value,and subject title (stitle)

**MUST DOWNLOAD**:
Please download the "Pipeline.py" into your working directory. There is no need to download the fastq files as examples because they will be downloaded through the script already. Any dependencies needed will also already automatically be downloaded through import statements at the top of the script. If they are not present, please use 

The example data download should already be included in the python script by using the wget and subprocess function. **If you need to change the data used go to line 51 and change the srr_numbers list and the urls in the download_fastq method**

"pip install "dependency-name" given the error thrown based on what's missing.

**How to run the script**:
Simply run the command python Pipeline.py -i /yourdirectory -o PipelineProject.log
