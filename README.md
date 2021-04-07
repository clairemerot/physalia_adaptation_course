# Physalia adaptation course

This repository includes scripts and data associated with the practical sessions of the Physalia Course about Adaptation Genomics
Prepared by Anna Tigano, Yann Dorant and Claire Mérot for the course in May 2021

Excpet for day 1 in which you won't have the bam files, all the tutorial can be completed with the compressed files stored in the github page for each day (day 2 to day 5


## Logging on the AWS server from your computer
### Mac OS X and Linux users
If you are using a Mac or Linux machine, you will need to open a terminal window and then use the ```ssh``` command. ssh stands for secure shell and is a way of connecting to and interacting with remote servers. You will need to log in to the cluster using both ssh and a keyfile that has been generated for you.
Firstly, download the keyfile and open a terminal window. Then copy the keyfile into your home directory like so:
```
cp mark.pem ~
```
Then, you should be able to log in with ssh. You need to provide ssh with the path to your key, which you can do with the ```-i``` flag. This basically points to your identity file or keyfile. For example:
```
ssh -i "~/anna.pem" anna@54.245.175.86
```
Note that you will need to change the log in credentials shown here (i.e. the username and keyfile name) with your own. Also be aware that the cluster IP address will change everyday. We will update you on this each day. You might be prompted to accept an RSA key - if so, just type yes and you will log in to the cluster!

### Downloading and uploading files
Occassionally, we will need to transfer files between the cluster and our local machines. To do this, we can use a command utility called ```scp```, which stans for secure copy. It works in a similar way to ssh. Let’s try making a dummy file in our local home directory and then uploading it to our home directory on the cluster.
```
# make a file
touch test_file
# upload to cluster
scp -i "~/anna.pem" test_file mark@54.245.175.86:~/
```
Just to break this down a little we are simply copying a file, ```test_file``` in this case to the cluster. After the `:` symbol, we are specifying where on the cluster we are placing the file, here we use `~/` to specify the home directory.


Copying files back on to our local machine is just as straightforward. You can do that like so:
```
# download to local
scp -i "~/mark.pem" mark@54.245.175.86:~/test_file ./
```
Where here all we did was use scp with the cluster address and path first and the location on our computers (in our working directory) second - i.e. `./`

Alternatively, we can use Cyberduck, or similar software, to transfer data back and forth.
Open Cyberduck, and click on 'Open connection' on the top left of the screen. Select 'SFTP (SSH File Tranfer Protocol)' from the dropdown menu and copy the IP address for the day in the Server field and your username for in the Username field. Leave the Password field empty, go down to SSH Private Key and add the Private Key Carlo sent you. Press Connect and you're connected!

### Windows users
[AWS for windows users](https://github.com/clairemerot/physalia_adaptation_course/blob/2021/AWS_windows.md)


## Day1: Handling NGS data: From raw reads to SNPs matrix

Capelin: Data taken from Cayuela et al,2020. Molecular Ecology https://doi.org/10.1111/mec.15499
 
Genome assembly: For this course, we made a dummy assembly of about 90 MB (instead of about 500 MB) and 5 chromosomes (instead of 24).

Raw reads: Obtained by GBS methods, sequenced with IonTorrent.
 
## Day2: Population structure and confounding factors

## Day3: Outlier detection - Environmental associations

## Day4: Accounting for Structural Variants

## Day5: Functional approaches

