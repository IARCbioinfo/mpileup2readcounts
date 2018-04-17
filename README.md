[![Build Status](https://travis-ci.org/gatoravi/mpileup2readcounts.svg?branch=master)](https://travis-ci.org/gatoravi/mpileup2readcounts)
## Synopsis
Get the readcounts at a locus by piping samtools mpileup output. The mpileup can contain one or several samples.
This program has been tested on samtools v1.3.1

## Install samtools 

## Compile mpileup2readcounts : 
```
g++ -std=c++11 -O3 mpileup2readcounts.cc -o mpileup2readcounts
```

## Usage

```
samtools mpileup -f ref.fa -l regions.bed BAM/*.bam | sed 's/		/	* 	*/g' | ./mpileup2readcounts 0 -5 false 3 0
```
Samtools arguments :
- FASTA file
- bed file
- BAM files : several samples can be parsed

Four options for mpileup2readcounts :
- 0 to parse all sample otherwise specify the number of the sample (for example 1 for the first sample)
- BQcut : base quality score cutoff for each mapped/unmapped base, only those larger than cutoff will be output in the result, to use no filter set BQcut to -5
- true to ignore indels
- min_ao : minimum number of non-ref reads in at least one sample to consider a site
- min_af : minimum allelic fraction in at least one sample to consider a site

## Example output


| chr |	loc	| ref	| depth	| A	| T	| C	| G	| a	| t	| c	| g	| Insertion	| Deletion	| depth	| A	| T	| C	| G	| a	| t	| c	| g	| Insertion	| Deletion	|
|-----------|--------------|-------------|-----------|--------------|-------------|-----------|--------------|-------------|-----------|--------------|-------------:|-----------|--------------|-------------|-----------|--------------|-------------|-----------|--------------|-------------|-----------|--------------|-------------|-----------|     
| 17 | 7572814	| C	| 28	| 0	| 3	| 23	| 0	| 0	| 0	| 2	| 0	| NA	| NA	| 8	| 0	| 0	| 8	| 0	| 0	| 0	| 0	| 0	| NA |	NA |
| 17	| 7572817	| C	| 32	| 2	| 0	| 26	| 0	| 2	| 0	| 2	| 0	| NA	| NA	| 8	| 0	| 0	| 8	| 0	| 0	| 0	| 0	| 0	| NA	| NA	|
| 17	| 7579643	| C	| 48	| 0	| 0	| 9	| 0	| 0	| 0	| 39	| 0	| NA	|	4:ccccagccctccaggt\|2:CCCCAGCCCTCCAGGT	| 9	| 0	| 0	| 6	| 0	| 0	| 0	| 3	| 0	| NA	| NA|


### Line content
Common information for all samples: 
- chromosome
- position on the chromosome
- reference base

For each sample :
- depth
- ATCG/atcg count
- insertions 
- deletions : in the example, for the first sample 6 deletions starting from position 7579643 + 1 are found

## Using GNU parallel tool
We first remove overlap in the bed using bedtools. Note that the header is removed here:
```
grep -v '^track' regions.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print $1"\t"$2"\t"$3}' | sed 's/[[:space:]]/:/' | sed 's/[[:space:]]/-/' | parallel --keep-order "samtools mpileup -f ref.fa --region {} test.bam | sed 's/		/	*	*/g' | mpileup2readcounts 0 -5 false 0 | tail -n +2"
```


