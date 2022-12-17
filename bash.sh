#!/bin/bash
# bash script for computational pipeline
# Inputs: test_id.txt ref_panel.txt iteration#

# Download Java for BEAGLE
module load java/1.8.0
# Create variable for BEAGLE
BEAGLE=/workdir/forensics/beagle/beagle.22Jul22.46e.jar
# Create variable for genetic map directory
MAP=/workdir/forensics/beagle/plink.GRCh37.map
# Create tmp directory for iteration
mkdir tmp${3}

# Create list of SNPsOnly vcf files
declare -a snpsonly
for file in /workdir/forensics/base_data/1kgp/vcf/*SNPsOnly.vcf
do
snpsonly=(${snpsonly[@]} $file)
done

# Create list of WithCODIS vcf files
declare -a withcodis
for file in /workdir/forensics/base_data/1kgp/vcf/*WithCODIS.vcf
do
withcodis=(${withcodis[@]} $file)
done

# Create list of CODIS loci
loci=(CSF1PO D13S317 D18S51 D3S1358 D5S818 D7S820 D8S1179 FGA TH01 TPOX vWA D1S1656 D2S441 D2S1338 D10S1248 D12S391 D19S433 D22S1045)

# Create associate array of CODIS loci to chromosome map file
declare -A chromosome
chromosome[CSF1PO]=/plink.chr5.GRCh37.map
chromosome[D13S317]=/plink.chr13.GRCh37.map
chromosome[D18S51]=/plink.chr18.GRCh37.map
chromosome[D3S1358]=/plink.chr3.GRCh37.map
chromosome[D5S818]=/plink.chr5.GRCh37.map
chromosome[D7S820]=/plink.chr7.GRCh37.map
chromosome[D8S1179]=/plink.chr8.GRCh37.map
chromosome[FGA]=/plink.chr4.GRCh37.map
chromosome[TH01]=/plink.chr11.GRCh37.map
chromosome[TPOX]=/plink.chr2.GRCh37.map
chromosome[vWA]=/plink.chr12.GRCh37.map
chromosome[D1S1656]=/plink.chr1.GRCh37.map
chromosome[D2S441]=/plink.chr2.GRCh37.map
chromosome[D2S1338]=/plink.chr2.GRCh37.map
chromosome[D10S1248]=/plink.chr10.GRCh37.map
chromosome[D12S391]=/plink.chr12.GRCh37.map
chromosome[D19S433]=/plink.chr19.GRCh37.map
chromosome[D22S1045]=/plink.chr22.GRCh37.map

# Step 1: Create vcf files of masked genotypes for individuals in the test set for each CODIS loci
for masked in ${snpsonly[@]}
do
tmp=${masked#*vcf/}
name=${tmp%_half*}
vcftools --vcf $masked --keep $1 --recode --out /workdir/forensics/tmp${3}/${name}maskedtest
done

# Step 2: Create vcf files of unmasked genotypes for individuals in the reference panel for each CODIS loci
for unmasked in ${withcodis[@]}
do
tmp=${unmasked#*vcf/}
name=${tmp%_half*}
vcftools --vcf $unmasked --keep $2 --recode --out /workdir/forensics/tmp${3}/${name}refdata
done

# Create vcf files of unmasked genotypes for individuals in the test set for each CODIS loci
for unmasked in ${withcodis[@]}
do
tmp=${unmasked#*vcf/}
name=${tmp%_half*}
vcftools --vcf $unmasked --keep $1 --recode --out /workdir/forensics/tmp${3}/${name}unmaskedtest
done

# Step 3: Extract GT of loci from unmasked genotypes for individuals in the test set for later comparison to imputed data
for str in ${loci[@]}
do
vcftools --vcf /workdir/forensics/tmp${3}/${str}unmaskedtest.recode.vcf --snp $str --extract-FORMAT-info GT --out /workdir/forensics/tmp${3}/${str}codistest
done

# Step 4: Use BEAGLE to phase the reference panel to generate SNP-STR haplotypes for each CODIS loci
for str in ${loci[@]}
do
java -Xmx1g -jar $BEAGLE gt=/workdir/forensics/tmp${3}/${str}refdata.recode.vcf out=/workdir/forensics/tmp${3}/${str}phased burnin=6 iterations=14 nthreads=1
done

# Step 5: Impute STR genotypes at CODIS loci of test set using BEAGLE and phased reference panel for each CODIS loci
for str in ${loci[@]}
do
chrom=${chromosome[$str]}
java -Xmx1g -jar $BEAGLE gt=/workdir/forensics/tmp${3}/${str}maskedtest.recode.vcf out=/workdir/forensics/tmp${3}/${str}imputed ref=/workdir/forensics/tmp${3}/${str}phased.vcf.gz gp=true impute=true map=${MAP}${chrom} nthreads=1
done

# Convert .vcf.gz files to .vcf files
gunzip tmp${3}/*.gz

# Create vcf file containing CODIS loci only of imputed test set (one file for each loci)
for str in ${loci[@]}
do
OUT=/workdir/forensics/tmp${3}/$str
vcftools --vcf /workdir/forensics/tmp${3}/${str}imputed.vcf --snp $str --recode --recode-INFO-all --out $OUT
done

# Step 6: Extract imputed genotype (GT) and genotype probabilities (GP) at each CODIS loci
for str in ${loci[@]}
do
FILE=/workdir/forensics/tmp${3}/${str}.recode.vcf
GTFILE=/workdir/forensics/tmp${3}/$str
GPFILE=/workdir/forensics/tmp${3}/$str
vcftools --vcf $FILE --extract-FORMAT-info GT --out $GTFILE
vcftools --vcf $FILE --extract-FORMAT-info GP --out $GPFILE
done
