capmq
=====

Cap the mapping quality in a given alignment file.

One way to deal with low level contamination in your alignment files is to cap
the mapping quality at some level that reflects your estimated certainty that
your reads do indeed come from your sample. If you esimate your level of
"contamination" (e.g. via a program like [VerifyBamID](http://genome.sph.umich.edu/wiki/VerifyBamID))
at, say `10^-3`, then your cannot say with more certainty than `MAPQ=30` that
your reads are mapped correctly:

	capmq -C30 in.bam

`capmq` uses [htslib](https://github.com/samtools/htslib) for SAM/BAM/CRAM
reading and writing.

```
git clone git://github.com/mcshane/capmq.git
cd capmq
make HTSDIR=/path/to/your/htslib/install/or/build/dir
make test
```
