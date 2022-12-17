# Cross-Database Matching of STR-SNP for Forensics

### Project Description
The objective of this project is to examine the effects of reference panel size on the imputation accuracy of CODIS loci for populations of different genetic ancestry. We used data from the 1000 Genomes Project that was unphased and split into a masked set that contained only SNPs and an unmasked set that contained SNPs and STRs. Next we used python to generate a random test set (of given population) and reference panel (of given size), and then used vcftools to extract the data for those individuals. Then, BEAGLE used this reference panel to phase SNP-STR haplotypes, and then impute the masked CODIS loci of the test set. We measured imputation accuracy by converting the genotypes to numpy arrays and computing the correlation coeffiecient for each codis loci. We tested imputation accuracy for five populations from various super-populations (Africans, Admixed Americans, East Asians, Europeans, and South Asians) and varied reference panel size from 10 to 50 individuals.

### Future Developments
Next steps would included running the pipeline with more populations, creating more complex reference panels, and running more iterations. Furthermore, we hope to improve the computational pipeline specifically for admixed populations. This will involve using a data set of simulated admixed genomes and running more imputation experiments to determine the optimum reference panel for admixed individuals.

### Installations
Running this pipeline requires downloading data from the 1000 genomes project. Additionally it requires the installation of BEAGLE ([beagle.22Jul22.46e.jar_](https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar)). Java 8 (https://www.java.com/en/download/manual.jsp) must be installed in order to use BEAGLE.

### System Requirements
Java 8
<br />BEAGLE
<br />vcftools
<br />Python
<br />Bash
<br />BioHPC

### References
Fortier AL, Kim J, Rosenberg NA, 2020. Human-Genetic Ancestry Inference and False Positives in Forensic Familial
Searching. G3 Genes|Genomes|Genetics, 10(8):2893–2902.

Saini S, Mitra I, Mousavi N, Fotsing SF, Gymrek M, 2018. A reference haplotype panel for genome-wide imputation of
short tandem repeats. Nature Communications, 9(1):4397.

Huang L, Li Y, Singleton AB, Hardy JA, Abecasis G, Rosenberg NA, Scheet P, 2009. Genotype-imputation accuracy across
worldwide human populations. The American Journal of Human Genetics, 84(2):235–250.

Edge MD, Algee-Hewitt BFB, Pemberton TJ, Li JZ, Rosenberg NA, 2017. Linkage disequilibrium matches forensic genetic
records to disjoint genomic marker sets. Proceedings of the National Academy of Sciences, 114(22):5671–5676.

Kim J, Edge MD, Algee-Hewitt BF, Li JZ, Rosenberg NA, 2018. Statistical detection of relatives typed with disjoint forensic
and biomedical loci. Cell, 175(3):848–858.e6.
4/5


