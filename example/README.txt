### Step 1 - For every sample, generate a strandcount.csv from the bam file ###
java -cp /your-path-to/HAROLD/lib/htsjdk-unspecified-SNAPSHOT.jar:/Your-path-to/HAROLD/lib/picocli-4.1.2.jar:/Your-path-to/HAROLD/lib/pal-1.5.1.jar:/Your-path-to/HAROLD/lib/cache2k-all-1.0.2.Final.jar:/Your-path-to/HAROLD/lib/commons-math3-3.6.1.jar:/Your-path-to/HAROLD/jar/MakeReadCount.jar makereadcount.MakeReadCount file.bam

### Step 2 - Running HaROLD ###
(This is an example command only, change the number of haplotypes ‘--haplotypes’ to get the optimised results)
java -jar /your-path-to/jar/Cluster_RG/dist/Cluster_RG.jar --count-file sample.txt --haplotypes 4 --alpha-frac 0.5 --gamma-cache 10000 -H -L --threads 4

### Step 3 - Refining output from HaROLD ###
(Input the result with the best likelihood from Step 2)
(Need to run this for every sample, here is just an example command for sample 2)
java -cp /your-path-to/lib/htsjdk-unspecified-SNAPSHOT.jar:/your-path-to/lib/picocli-4.1.2.jar:/your-path-to/lib/pal-1.5.1.jar:/your-path-to/lib/commons-math3-3.6.1.jar:/your-path-to/lib/cache2k-all-1.0.2.Final.jar:/your-path-to/lib/flanagan.jar:/your-path-to/jar/RefineHaplotypes.jar refineHaplotypes.RefineHaplotypes -t sample2 --bam sample2-1longitudinal.sorted.dedup.bam.fixed.bam --baseFreq nhaplo_4_results.lld --refSequence refseq-JX459907.fasta --hapAlignment nhaplo_4_resultsHaplo.fasta --iterate

