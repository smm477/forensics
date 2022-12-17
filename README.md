# Cross-Database Matching of STR-SNP for Forensics

### Project Description
The objective of this project is to examine the effects of reference panel size on the imputation accuracy of CODIS loci for populations of different genetic ancestry. We used data from the 1000 Genomes Project that was unphased and split into a masked set that contained only SNPs and an unmasked set that contained SNPs and STRs. Next we used python to generate a random test set (of given population) and reference panel (of given size), and then used vcftools to extract the data for those individuals. Then, BEAGLE used the reference panel to phase SNP-STR haplotypes, and then impute the masked CODIS loci of the test set. Then, we measured imputation accuracy by converting the genotypes to numpy arrays and computing the correlation coeffiecient for each codis loci. We tested imputation accuracy for five populations from various super-populations (Africans, Admixed Americans, East Asians, Europeans, and South Asians) and varied reference panel size from 10 to 50 individuals.

### Future Developments
Next steps would included running the pipeline with more populations, creating more complex reference panels, and running more iterations. Furthermore, we hope to improve the computational pipeline specifically for admixed populations. This will involve using a data set of simulated admixed genomes and running more imputation experiments to determine the optimum reference panel for admixed individuals.

### Software Used
Java 8
<br />BEAGLE
<br />vcftools
<br />Python
<br />BioHPC

