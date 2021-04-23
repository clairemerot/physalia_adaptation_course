Open R in your terminal and let's simulate genotypes
```
args = commandArgs(trailingOnly=TRUE)

#source functions from Gautier file
source("/home/ubuntu/00-softwares/baypass_2.2/utils/baypass_utils.R")

#load omega matrix
omega = as.matrix(read.table("05_baypass/prunedsnps.output_mat_omega.out"))

#load beta parameters
pi.beta.coef=read.table("05_baypass/allsnps.output_summary_beta_params.out", h=T)$Mean

#create the POD - simulate
simu.bta <- simulate.baypass(omega.mat=omega, nsnp=500,beta.pi=pi.beta.coef, suffix="simulates_pods")
```
You can now quit R, mv the simulated file into the folder baypass and run baypass on the simulated data 
```
mv G.simulates_pods 05_baypass/
```
