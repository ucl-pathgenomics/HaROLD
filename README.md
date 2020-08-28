# HaROLD - HAplotype Reconstruction Of Longitudinal Deep sequencing data

## Overview
HaROLD reconstructs haplotypes based on identifying co-varying variant frequencies using a probabilistic framework. For more details, please refer to [our preprint](https://www.biorxiv.org/content/10.1101/444877v2) on bioRxiv.

## Usage

#### Step 1 - BAM Files Conversion
Convert the BAM files using [Samtools](http://www.htslib.org).
```sh
samtools view -h -G69 to.convert.bam | samtools view -h -G133 > file.bam
```
For every sample, generate a strandcount.csv from the BAM file.
```sh
java -cp /your-path-to/HAROLD/lib/htsjdk-unspecified-SNAPSHOT.jar:
/Your-path-to/HAROLD/lib/picocli-4.1.2.jar:
/Your-path-to/HAROLD/lib/pal-1.5.1.jar:
/Your-path-to/HAROLD/lib/cache2k-all-1.0.2.Final.jar:
/Your-path-to/HAROLD/lib/commons-math3-3.6.1.jar:
/Your-path-to/HAROLD/jar/MakeReadCount.jar makereadcount.MakeReadCount file.bam
```

#### Step 2 - Running HaROLD
Change the number of haplotypes ‘--haplotypes’ to get the optimised results.
```sh
java -jar /your-path-to/jar/Cluster_RG/dist/Cluster_RG.jar
--count-file sample.txt --haplotypes 4 --alpha-frac 0.5 --gamma-cache 10000 -H -L --threads 4
```

#### Step 3 - Refining Output from HaROLD
Input result with the best likelihood from Step 2.
This is an example for sample 2. Run this for every sample.
```sh
java -cp /your-path-to/lib/htsjdk-unspecified-SNAPSHOT.jar:
/your-path-to/lib/picocli-4.1.2.jar:
/your-path-to/lib/pal-1.5.1.jar:
/your-path-to/lib/commons-math3-3.6.1.jar:
/your-path-to/lib/cache2k-all-1.0.2.Final.jar:
/your-path-to/lib/flanagan.jar:
/your-path-to/jar/RefineHaplotypes.jar refineHaplotypes.RefineHaplotypes
-t sample2 --bam sample2-1longitudinal.sorted.dedup.bam.fixed.bam
--baseFreq nhaplo_4_results.lld --refSequence refseq-JX459907.fasta
--hapAlignment nhaplo_4_resultsHaplo.fasta --iterate
```

## Getting help
If you have any questions, please feel free to contact [Juanita Pang](mailto:juanita.pang.16@ucl.ac.uk) or [Richard Goldstein](mailto:r.goldstein@ucl.ac.uk).
