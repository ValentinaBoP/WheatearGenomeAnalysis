# De-novo characterisation of repetitive elements

## RepeatModeler2 + classifier (TBD)

### RepeatModeler2

```
# go to working directory and create folder for the de novo characterisation
cd /proj/sllstore2017073/private/Valentina/2020Wheatear
mkdir -p Intermediate/DenovoTE
cd Intermediate/DenovoTE
```

Link reference assembly

```
ln -s /proj/sllstore2017073/private/Octavio/Wheatear_assembly/3DDNA_new_curated_assembly_Pilon_2_final_reference/OenMel_1.1_complete.fasta ./oenMel.fasta
```

Run RepeatModeler2 on `oenMel.fasta`

```
cd /proj/sllstore2017073/private/Valentina/2020Wheatear/Code
```

Save the following SLURM job into `Code` as `RMDL_oenMel.sh`

```
#!/bin/bash -l
#SBATCH -J RMDL_oenMel
#SBATCH -o RMDL_oenMel.output
#SBATCH -e RMDL_oenMel.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 72:00:00
#SBATCH -A snic2021-5-585
#SBATCH -n 16

# load modules
ml bioinfo-tools RepeatModeler/2.0.1

# go to working directory
cd /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/DenovoTE

# run RepeatModeler2
# STEP1: create a database for the genome
BuildDatabase -name oenMel oenMel.fasta

# STEP2: run RepeatModeler2 on the database you just created
RepeatModeler -database oenMel -LTRStruct -pa 4
```

Submit job

```
sbatch RMDL_oenMel.sh
```

Remove the intermediate files from RepeatModeler2 (it may take a while)

```
cd /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/DenovoTE
rm -r RM_23941.TueMar221108112022
```

### Rename the consensus sequences
Rename the raw consensus sequences with the species abbreviation 'oenMel'

```
cd /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/DenovoTE
perl /proj/sllstore2017073/private/scripts/renameRMDLconsensi_Vale.pl oenMel-families.fa oenMel oenMel_rm1.0.fasta
```

### Remove protein coding genes from oenMel_rm1.0
Input: `oenMel_rm1.0.fasta`
Output: `oenMel_rm1.1.fasta`

```
cd /proj/sllstore2017073/private/Valentina/2020Wheatear/Tools/repeatlib_filtering
tmux new-session -s lib_filter
conda deactivate
conda activate sm5191py36
snakemake -j 100 --use-conda --cluster-config cluster.yaml --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time}" &> main_run.out


snakemake -j 100 --cluster-config cluster.yaml --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time}" &> main_run.out
# detach from tmux
```

### Run RepeatMasker with oenMel_rm1.0

```
#!/bin/bash -l
#SBATCH -J RMSK_oenMel_rm1.0
#SBATCH -o RMSK_oenMel_rm1.0.output
#SBATCH -e RMSK_oenMel_rm1.0.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 72:00:00
#SBATCH -A snic2021-5-585
#SBATCH -n 20

# load modules
ml bioinfo-tools RepeatMasker/4.1.0

# go to working directory
mkdir /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/oenMel_rm1.0
OUTDIR=/proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/oenMel_rm1.0

# copy files into temporary directory
cp /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/DenovoTE/oenMel_rm1.0.fasta $SNIC_TMP
cp /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/DenovoTE/oenMel.fasta $SNIC_TMP
cd $SNIC_TMP

# run RepeatMasker
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -s -dir $OUTDIR -lib oenMel_rm1.0.fasta oenMel.fasta
```

### Create divsum file and CpG-corrected out file

```
#!/bin/bash -l
#SBATCH -J RMSK_oenMel_rm1.0
#SBATCH -o RMSK_oenMel_rm1.0.output
#SBATCH -e RMSK_oenMel_rm1.0.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 06:00:00
#SBATCH -A snic2021-5-585
#SBATCH -n 2

# go to working directory
cd /home/vpeona/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/oenMel_rm1.0

ALIGN=oenMel.fasta.align
DIVSUM=${ALIGN%.*}.divsum

# Convert the .align file into a table with Kimura 2-parameter distances (excluding CpG sites)
perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG $ALIGN > $ALIGN.k2p.noCpG

# Add sizes of the individual TE copies
perl /proj/sllstore2017073/private/scripts/RM_divergence/RM.addsize.div.pl $ALIGN.k2p.noCpG > $ALIGN.k2p.noCpG.size

# Remove the term "kimura=" from the last column
sed -i 's/kimura=//g' $ALIGN.k2p.noCpG.size

# load modules
ml bioinfo-tools RepeatMasker/4.1.2-p1

calcDivergenceFromAlign.pl -s $DIVSUM -noCpGMod $ALIGN

createRepeatLandscape.pl -div $DIVSUM
```

### Run RepeatMasker with oenMel_rm1.1

```
#!/bin/bash -l
#SBATCH -J RMSK_oenMel_rm1.1
#SBATCH -o RMSK_oenMel_rm1.1.output
#SBATCH -e RMSK_oenMel_rm1.1.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 72:00:00
#SBATCH -A snic2021-5-585
#SBATCH -n 20

# load modules
ml bioinfo-tools RepeatMasker/4.1.0

# go to working directory
mkdir /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/oenMel_rm1.1
OUTDIR=/proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/oenMel_rm1.1

# copy files into temporary directory
cp /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/DenovoTE/oenMel_rm1.1.fasta $SNIC_TMP
cp /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/DenovoTE/oenMel.fasta $SNIC_TMP
cd $SNIC_TMP

# run RepeatMasker
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -s -dir $OUTDIR -lib oenMel_rm1.1.fasta oenMel.fasta
```

### Run RepeatMasker with oenMel_rm1.1 + avian TE library

```
#!/bin/bash -l
#SBATCH -J RMSK_avianLib
#SBATCH -o RMSK_avianLib.output
#SBATCH -e RMSK_avianLib.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 72:00:00
#SBATCH -A snic2021-5-585
#SBATCH -n 20

# load modules
ml bioinfo-tools RepeatMasker/4.1.0

# go to working directory
mkdir /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/avianTE
OUTDIR=/proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/avianTE

# copy files into temporary directory
cp /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/DenovoTE/bird_library_25Oct2020_oenMel_rm1.1.fasta $SNIC_TMP
cp /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/DenovoTE/oenMel.fasta $SNIC_TMP
cd $SNIC_TMP

# run RepeatMasker
RepeatMasker -pa 20 -a -xsmall -gccalc -excln -s -dir $OUTDIR -lib bird_library_25Oct2020_oenMel_rm1.1.fasta oenMel.fasta
```


### Create divsum file and CpG-corrected out file

```
#!/bin/bash -l
#SBATCH -J CALCDIV_oenMel_rm1.1
#SBATCH -o CALCDIV_oenMel_rm1.1.output
#SBATCH -e CALCDIV_oenMel_rm1.1.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 06:00:00
#SBATCH -A snic2021-5-585
#SBATCH -n 2

# go to working directory
cd /home/vpeona/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/oenMel_rm1.1

ALIGN=oenMel.fasta.align
DIVSUM=${ALIGN%.*}.divsum

# Convert the .align file into a table with Kimura 2-parameter distances (excluding CpG sites)
perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG $ALIGN > $ALIGN.k2p.noCpG

# Add sizes of the individual TE copies
perl /proj/sllstore2017073/private/scripts/RM_divergence/RM.addsize.div.pl $ALIGN.k2p.noCpG > $ALIGN.k2p.noCpG.size

# Remove the term "kimura=" from the last column
sed -i 's/kimura=//g' $ALIGN.k2p.noCpG.size

# load modules
ml bioinfo-tools RepeatMasker/4.1.2-p1

calcDivergenceFromAlign.pl -s $DIVSUM -noCpGMod $ALIGN

createRepeatLandscape.pl -div $DIVSUM
```

```
#!/bin/bash -l
#SBATCH -J CALCDIV_avianLib
#SBATCH -o CALCDIV_avianLib.output
#SBATCH -e CALCDIV_avianLib.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 06:00:00
#SBATCH -A snic2021-5-585
#SBATCH -n 2

# go to working directory
cd /home/vpeona/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/avianTE

ALIGN=oenMel.fasta.align
DIVSUM=${ALIGN%.*}.divsum

# Convert the .align file into a table with Kimura 2-parameter distances (excluding CpG sites)
perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG $ALIGN > $ALIGN.k2p.noCpG

# Add sizes of the individual TE copies
perl /proj/sllstore2017073/private/scripts/RM_divergence/RM.addsize.div.pl $ALIGN.k2p.noCpG > $ALIGN.k2p.noCpG.size

# Remove the term "kimura=" from the last column
sed -i 's/kimura=//g' $ALIGN.k2p.noCpG.size

# load modules
ml bioinfo-tools RepeatMasker/4.1.2-p1

calcDivergenceFromAlign.pl -s $DIVSUM -noCpGMod $ALIGN

createRepeatLandscape.pl -div $DIVSUM > $DIVSUM.html
```

## Repeat landscape

```{bash}
# RMSK annotation file
RMSK=oenMel.fasta.align.k2p.noCpG.size

# Index file
INDEX=../../../Data/oenMel.fasta.fai

# load R libraries
ml R_packages/4.1.1

# run script to get the general landscape for the entire genome
# arguments of the script: mode, annotation file, name of the plot file
# the mode genome plots one general landscape for the entire genome (all sequences)
Rscript --vanilla $CODE/landscapes1.2.R genome $RMSK oenMel_tot_landscape

# run script to get the general landscape for the entire genome
# arguments of the script: mode, annotation file, name of the plot file
# the mode genome plots one general landscape for the entire genome (all sequences)
Rscript --vanilla $CODE/landscapes1.2.R genome percentage $RMSK $INDEX oenMel_tot_landscape_percentage
```

## Repeat pie-chart

I use the .tbl file to build the pie-chart

```{bash}
cd /proj/sllstore2017073/private/Valentina/2020Wheatear/Intermediate/RMSK/avianTE

printf "SINE\t0.16\nLINE\t4.01\nLTR\t4\nDNA\t0.32\nUnknown\t1.19\nNon-repetitive DNA\t90.32\n" > pieChart_data

Rscript --vanilla $CODE/pieChart1.0.R
```










