# Day1: Handling NGS data: From raw reads to SNPs matrix

## Step 1 : XX

### dacdca

Welcome to the practical portion of the Adaptation Genomics course!

Here are the instructions to log in the AWS server.

Once you have successfully logged in, the first thing you will do is to create your own folder and generate a symbolic link so that we don't need to personalize all the scripts
```
mkdir *yourlastname*
ln -s ~/*yourlastname* ~/home
```

Copy this repository in your folder
```
cd ~/home
git clone  https://github.com/clairemerot/physalia_adaptation_course.git
```


## Step 1 : Mapping short reads to the capelin (dummy) reference genome
The data we'll be working on are from capelin (Mallotus villosus), a small marine fish, from Cayuela et al. 2020, Molecular Ecology. They were generated on the IonTorrent platform but all the preprocessing (demultiplexing, adapter removal, read trimming, quality filtering) has already been done. You can find example script for these preliminary steps here and here.

To run analysis that take some time we will submit jobs through bash scripts.
To map the filtered data to the capelin reference genome, make the script ```bowtie_mapping.sh``` in folder xxx executable and submit it with
```
chmod +x bowtie_mapping.sh
./bowtie_mapping.sh
```
To run a bash script like this you can also type
```
bash bowtie_mapping.sh
``` 
without making it executable.

