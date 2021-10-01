# HaROLD - HAplotype Reconstruction Of Longitudinal Deep sequencing data
<p align="center">
<img src="https://github.com/ucl-pathgenomics/HaROLD/blob/master/HaROLD_logo.png" width="350">

## Overview
HaROLD reconstructs haplotypes based on identifying co-varying variant frequencies using a probabilistic framework. For more details, please refer to [our preprint](https://www.biorxiv.org/content/10.1101/444877v2) on bioRxiv.

## Usage

#### Prepare input files
This step might not be necessary depending on software used for alignment. Tested with bam files produced by Bbmap and bwa. 
Convert the BAM files using [Samtools](http://www.htslib.org).

```sh
samtools view -h -G69 to.convert.bam | samtools view -h -G133 > file.bam
```

For every sample, generate a strandcount.csv from the BAM file.

```sh
java -cp /your-path-to/HAROLD/lib/htsjdk-unspecified-SNAPSHOT.jar: \
/Your-path-to/HAROLD/lib/picocli-4.1.2.jar: \
/Your-path-to/HAROLD/lib/pal-1.5.1.jar: \
/Your-path-to/HAROLD/lib/cache2k-all-1.0.2.Final.jar: \
/Your-path-to/HAROLD/lib/commons-math3-3.6.1.jar: \
/Your-path-to/HAROLD/jar/MakeReadCount.jar makereadcount.MakeReadCount file.bam
```

#### Step 1 - Running HaROLD
HaROLD requires as input file only a list of strandcount.csv files. Longitudinal data are listed together, but separate samples from different runs need to be submitted separately. 
For example, if you have 4 longitudinal samples from the same patient, your sample.txt will look like:

```sh
sampleA_timepoint1.strandcount.csv
sampleA_timepoint2.strandcount.csv
sampleA_timepoint3.strandcount.csv
sampleA_timepoint4.strandcount.csv
```

HaROLD can be run as: 

```sh
java -jar /your-path-to/HaROLD/jar/Cluster_RG/dist/Cluster_RG.jar --help
```

Available parameters for HaROLD can be found in : 

```sh

Usage: richards-haplotype-model [-AhHLNvV] [--alpha-frac=<alpha_frac>]
                                [--threads=<threads>] [--tol=<tol>]
                                [-o=<optimiser>] [-p=<prefix>]
                                [-s=<randomSeed>] [-a=<initialAlphaParams>
                                <initialAlphaParams>]... -c=<countFile>...
                                [-c=<countFile>...]...
                                [-f=<initialFreqFile>...]... -n=<haplotypes>...
                                [-n=<haplotypes>...]...

      --alpha-frac=<alpha_frac>
                            Fraction of sites to use to optimise error parameters
      --threads=<threads>
      --tol=<tol>
  -a, --initial-alpha=<initialAlphaParams> <initialAlphaParams>

  -A, --fix-alpha           Fix alpha parameters
  -c, --count-file=<countFile>...
                            file containing list of count files
  -f, --initial-freq-file=<initialFreqFile>...
                            optional file containing hap frequency values
  -h, -?, --help            give this help list
  -H, --printHaplotypes     Print haplotypes
  -L, --printLikelihoods    Print likelihoods
  -n, --haplotypes=<haplotypes>...
                            number of haplotypes
  -N, --noOpt               Process without optimising
  -o, --optimiser=<optimiser>
                            Optimiser for haplotype frequencies
  -p, --prefix=<prefix>     Results file prefix
  -s, --seed=<randomSeed>
  -v, --verbose
  -V, --version
```

For example, the following command:


```sh
java -jar /your-path-to/jar/Cluster_RG/dist/Cluster_RG.jar \
--count-file sample.txt --haplotypes 4 --alpha-frac 0.5 --gamma-cache 10000 \
-H -L --threads 4 -p /your-path-to-results/Step1_results
```

The output of this step: 

- your-path-to-results/Step1_results.log: log file 
- your-path-to-results/Step1_results.fasta: haplotypes fasta sequences 
- your-path-to-results/Step1_results.lld: base frequency file 

NB. How to choose the number of haplotypes: the log files provides the total likelihood for the model (at the end of log file, "Main: Final total likelihood"), you should choose the number of haplotypes that gives you the highest likehood. However, the following step will help to find the best number. 

#### Step 2 - Refining Output from HaROLD
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
