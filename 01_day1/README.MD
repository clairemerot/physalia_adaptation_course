# Day1: Handling NGS data: From raw reads to SNPs matrix

Welcome to the practical portion of the Adaptation Genomics course!

## Step 1 : Logging on the AWS server from your computer
### Mac OS X and Linux users
If you are using a Mac or Linux machine, you will need to open a terminal window and then type ```ssh```. ssh stands for secure shell and is a way of interacting with remote servers. You will need to log in to the cluster using both ssh and a keyfile that has been generated for you.
Firstly, download the keyfile and open a terminal window. Then copy it into your home directory like so:
```
cp mark.pem ~
```
Then you should be able to log in with ssh whatever your working directory is. You need to provide ssh with the path to your key, which you can do with the ```-i``` flag. This basically points to your identity file or keyfile. For example:
```
ssh -i "~/anna.pem" anna@54.245.175.86
```
Note that you will need to change the log in credentials shown here (i.e. the username and keyfile name) with your own. Also be aware that the cluster IP address will change everyday. We will update you on this each day. You might be prompted to accept an RSA key - if so just type yes and you will log in to the cluster!

###Downloading and uploading files
Occassionally, we will need to shift files between the cluster and our local machines. To do this, we can use a command utility called ```scp``` or secure copy. It works in a similar way to ssh. Let’s try making a dummyfile in our local home directory and then uploading it to our home directory on the cluster.
```
# make a file
touch test_file

# upload to cluster
scp -i "~/anna.pem" test_file mark@54.245.175.86:~/
```
Just to break this down a little we are simply copying a file, ```test_file``` in this case to the cluster. After the : symbol, we are specifying where on the cluster we are placing the file, here we use ~/ to specify the home directory.

Copying files back on to our local machine is just as straightforward. You can do that like so:
```
# download to local
scp -i "~/mark.pem" mark@54.245.175.86:~/test_file ./
```
Where here all we did was use scp with cluster address first and the the location (our working directory) second - i.e. ./

Alternatively, we can use Cyberduck, or similar software, to transfer data back and forth.
Open Cyberduck, and click on 'Open connection' on the top left of the screen. Select 'SFTP (SSH File Tranfer Protocol)' from the dropdown menu and copy the IP address for the day in the Server field and your username for in the Username field. Leave the PAssword field empty, go down to SSH Private Key and add the Private Key Carlo sent you. Press Connect and you're connected!

### Windows users
If you are using a Windows machine, you will need to log on using PuTTY since there is no native ssh client. PuTTY does not natively support the private key format (.pem) needed to login to our Amazon cloud instance. You first need to convert the private key that we gave to you to a key that PuTTY can read. PuTTY has a tool named PuTTYgen, which can convert keys to the required PuTTY format (.ppk). When you installed PuTTY, it will also have installed PuTTYgen.

First, start PuTTYgen (for example, from the Start menu, choose All Programs > PuTTY > PuTTYgen). Then select RSA and click on Load:

In the new window that pops up, Change “PuTTY Private Key Files” to “All Files” to allow you to find your pem file.

Then save your key and click on YES to dismiss the Warning as shown below.

Great, now your key file is ready and we can start Putty. In Putty, enter your user name and the IP address in the format <user_name>@<IP adress>. Make sure that 22 is given as Port and that SSH is selected.

Next, on the menu on the left, expand the “SSH” panel by clicking on the + and select “Auth”. Then, select your new putty format key file (.ppk) with Browse. Do NOT click on Open yet.

To make sure that you will not get logged out if you are inactive for a while, you can configure PuTTY to automatically send ‘keepalive’ signals at regular intervals to keep the session active. Click on Connection, then insert 180 to send keepalive signals every 3 minutes. Do NOT click on Open yet.

To avoid having to change the settings each time you log in, you can save the PuTTY session. Click on Session to get back to the basic options window. Once you are happy with all settings, give the session a name and click Save. Now you are ready to start the session with “Open”. The first time PuTTY will display a security alert dialog box that asks whether you trust the host you are connecting to. Click yes.

When you log in the next time, you can just click on the saved session and click load. If the IP address changed in the mean time (e.g. because we stopped the Amazon instance over night), you will need to replace the IP address by the new one. I would then recommend to Save the changed settings. Then simply click Open to start the session.

If the IP address did not change and you just want to login again, you can also right-click on the putty symbol in the taskbar (provided that you have pinned it to the taskbar) and select the session.

You can tranfer files back and forth with Filezilla, a handy software to move files from a remote server such as the Amazon cloud or a cluster of your university.

Open Filezilla and choose Edit -> Settings.
Next, choose SFTP and Add the .pem key file as indicated below and click OK.
Finally, enter the IP address and the user name and when you hit enter, it should connect you. Next, time, you can use the Quickconnect dropdown menu, provided the IP address has not changed in the meantime.
Now you will see the file directory system (folders) on your local computer on the left and your folders on the amazon cloud on the right. You can now just drag and drop files from one side to the other.

## Step 2: Getting familiar with the UNIX environment
Although you should have familiarized with working on a UNIX environment before the beginning of the course, here below are the commands that we will use most commonly. After you have logged in, create your own folder and move in it. You will run all your scripts from here
```
mkdir yourname
cd yourname
```
Then copy some of the files we will need today in a new folder
```
mkdir scripts
cd scripts
cp ~/Share/physalia_adaptation_course/00_documents/* ./ ###this command copies all the files in the folder 00_documents into the folder you're currently in
```
Use ```ls``` to get a list of the files we have copied over, and look at the content of one of those files with
```
less popmap_2lin.txt 
```
You can stroll down by pressing the bar on your keyboard and ```q``` to exit from the screen.

If you want to print a few lines to the screen (which won't disapper) type
```
head popmap_2lin.txt
```
which will print the first 10 lines by default. You can specify the number of lines with ```-n 25```. Same for ```tail```, which shows the end of the file.

To edit scripts and files, we'll use ```nano``` as text editor
Type 
```
nano popmap_2lin.txt
``` 
and move around using the arrows on your keyboard.
At the bottom of your screen are the shortcuts for different commands. For example, to modify the file and save it with a different name, write something in the first line, press ctlr + x and edit the name of the file into ```test_nano.txt```. Hit enter and the file will be saved.
Now, let's remove this file with
```
rm test_nano.txt
```
WARNING!!! ```rm``` removes files and folders for good, and they can't be retrieved. So please keep that in mind!

## Step 3: Getting from raw data to mapped data
The data we'll be working on are from capelin (Mallotus villosus), a small marine fish, from Cayuela et al. 2020, Molecular Ecology. They were generated on the IonTorrent platform but all the preprocessing (demultiplexing, adapter removal, read trimming, quality filtering) has already been done. You can find example script for these preliminary steps in the scripts folder of day 1.

As we're analyzing a large dataset that would take many hours to align to the reference genome, we also did the mapping for you using this script ```physalia_adaptation_course/01_day1/01_scripts/utilities/04_bwa_mem_align_ionproton.sh```. 
Here I will provide you with a step-by-step breakdown of that script.
The file starts with the line
```
#!/bin/bash
```
We will write your scripts in bash, and that is the line that tell the server in what language the script is written in.
A few lines starting with ```#SBATCH``` follow. We won't need these lines because we won't use a job scheduler, but this is a good example for the people that have SLURM on their institution's servers.

We set a few variables as shown below, These are very handy to make this script easy to follow and more customizable. For example, the number of CPUs is not hard-coded in the script and can be conveniently changed when submitting the script.
```
# Global variables
GENOMEFOLDER="02_genome"
GENOME="genome_mallotus_dummy.fasta"
DATAFOLDER="03_raw_reads"
ALIGNEDFOLDER="04_aligned_files"
NCPU=$4
```

First, it uses a for loop in bash to loop through the list of fastq files to map to the reference genome
```
for file in $(ls -1 "$DATAFOLDER"/*.fq.gz)
do
```
then, it uses a one-liner (even though it looks like multiple lines for the use of backslashes \) to map the data and pipe (with |) the .sam output into the binary format .bam
```
bwa mem -t "$NCPU" -k 19 -c 500 -O 0,0 -E 2,2 -T 0 \
        -R "$ID" \
        "$GENOMEFOLDER"/"$GENOME" "$DATAFOLDER"/"$name" 2> /dev/null |
        samtools view -Sb -q 1 -F 4 -F 256 -F 2048 \
        - > "$DATAFOLDER"/"${name%.fq.gz}".bam
```
The alignments in bam files are sorted and indexed
```
samtools sort --threads "$NCPU" -o "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam \
        "$DATAFOLDER"/"${name%.fq.gz}".bam

    samtools index "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam
```
and the intermediate files are deleted with
```
rm "$DATAFOLDER"/"${name%.fq.gz}".bam
```



mkdir scripts
--> here go scripts and popmaps
cp ~/Share/physalia_adaptation_course/01_day1/01_scripts/\*.sh scripts/
cp ~/Share/physalia_adaptation_course/00_documents/* scripts/


## Step 1 : Mapping short reads to the capelin (dummy) reference genome
The data we'll be working on are from capelin (Mallotus villosus), a small marine fish, from Cayuela et al. 2020, Molecular Ecology. They were generated on the IonTorrent platform but all the preprocessing (demultiplexing, adapter removal, read trimming, quality filtering) has already been done. You can find example script for these preliminary steps here and here.

To run analysis that take some time we will submit jobs through bash scripts. Below is an example of a short bash scripts
```
#!/bin/bash


setsid nohup serve.js &
nohup ./model >outfile.txt &

Mapping is our computational bottleneck we have done it for you using this script here https://github.com/enormandeau/stacks_workflow/blob/master/00-scripts/bwa_mem_align_reads.sh


        
To run analysis that take some time we will submit jobs through bash scripts.
To map the filtered data to the capelin reference genome, make the script ```bowtie_mapping.sh``` in folder xxx executable and submit it with
```
chmod +x physalia_adaptation_course/01_day1/01_scripts/bowtie_mapping.sh
physalia_adaptation_course/01_day1/01_scripts/bowtie_mapping.sh physalia_adaptation_course/00_documents/samplelist_all.txt
```
To run a bash script like this you can also type
```
bash physalia_adaptation_course/01_day1/01_scripts/bowtie_mapping.sh physalia_adaptation_course/00_documents/samplelist_all.txt
``` 
without making it executable.

### File formats
While we wait we can explore the read file (.fastq format) and the alignment file (.sam/.bam format), which are the input and output of the alignemnt, respectively.
#### FASTQ
```
@NB551191:35:HNWGCBGX3:1:11101:9622:1037 1:N:0:AGGCAGAA+ACTCTAGG
CCTTGNTGCACCTGTGGCATGGAGACCGAATCTTGTGGGGAAACAATCATTTCTTCAGGTCTGAGCTCTCAGATTT
+
AAAAA#EAEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEAE/E/AEEEEEEEEEE/EEEEEEEEEEEEEEEEEEE
@NB551191:35:HNWGCBGX3:1:11101:20086:1038 1:N:0:AGGCAGAA+ACTCTAGG
CTGGTNCACCATCCTTGTGTGCTGTTTCATGACAGTAATTACTGAGAGGGTCTGCAATTCAGATCACCTGAAACTC
+
AAAAA#EEEEAEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@NB551191:35:HNWGCBGX3:1:11101:4269:1051 1:N:0:AGGCAGAA+ACTGTAGG
TAACAGCACAGAGGATTGAATAAGGTGAGAGCAAAAGTCCTACTACTTATTCAGGCCCCATGTAGCAGTATTCCTC
+
AAAAA6EEEEEEEEEEEEEA/EEEEEEEAEEEEEEE6/EA/EEEEEE/E//EA//EEEEEEEE6A/E6/EE/EEEE
```
The snippet above shows 3 reads, with each read coming with 4 lines of information.
1. Sequence ID and information
2. The actual read sequence
3. Generally plus a '+' sign but it could contain additional information
4. The quality scores for each base in line 2

Knowing this, we can calculate the number of reads in a file by counting the number of lines
```
wc -l file.fastq
```
and divide by 4.

When we have paired reads, they are generally stored in two separate files and are paired by position: the first read in ```read1.fastq``` is paired with the first read in ```read2.fastq```, the second read in ```read1.fastq``` is paired with the second read in ```read2.fastq```, and so on. One way to check if your files are intact is to ensure that the two read files have the same amount of lines/reads.
#### SAM/BAM
The alignment file contains much more information. In addition to all the information contained in the .fastq files, we now have information of the quality of the alignment and the position of where those reads mapped in the genome. Although we produce .sam files with the alignment, we quickly convert them to .bam, which is in binary format. These files are not readable directly (but you can use ```samtools view```).

Back to mapping...

