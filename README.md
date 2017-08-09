PRINIA
=======

Prinia is a python package for designing primers from genomic regions or LOVD import format files. 

Requirements
-------------

You must have a recent version of [primer3](http://primer3.sourceforge.net/), [samtools](http://www.htslib.org/download/) and [bwa](https://github.com/lh3/bwa)

Furthermore, the following python packages are required:

* future >= 0.15.2
* biopython >= 1.66
* lxml >= 3.4.4
* pyfaidx >= 0.4.4
* fastools >= 0.12.0
* pyvcf >= 0.6.7 


Installation
-------------

It is recommended you use a [virtual environment](https://virtualenv.readthedocs.org), such that the risk of dependency hell is minimized.
  
After creating your virtualenv, clone this repository. 
Install all the required python packages with:
```
pip install -r requirements.txt
```


Then install PRINIA with:

```
python setup.py develop
```


Primer design
-------------

A `primerdesign` tool will be be added to your path upon installation of prinia.
This tool will generate your primers from BED records (regions) or LOVD import file format files.   
 
### Help


```
usage: primerdesign [-h] (-l LOVD | --region REGION) [-p PADDING]
                    (-x XML | -t TSV) -b BAM [-s SAMPLE]
                    [--product_size PRODUCT_SIZE] [--min-margin MIN_MARGIN]
                    [--strict] [--n_raw_primers N_RAW_PRIMERS] [--m13]
                    [--m13-forward M13_FORWARD] [--m13-reverse M13_REVERSE]
                    [-f FIELD] [-af ALLELE_FREQ] [-fq1 FQ1] [-fq2 FQ2] -R
                    REFERENCE --dbsnp DBSNP --primer3 PRIMER3 [--bwa BWA]
                    [--samtools SAMTOOLS] [--ignore-errors]
                    [--opt-primer-length OPT_PRIMER_LENGTH]
                    [--opt-gc-perc OPT_GC_PERC]
                    [--min-melting-temperature MIN_MELTING_TEMPERATURE]
                    [--max-melting-temperature MAX_MELTING_TEMPERATURE]

optional arguments:
  -h, --help            show this help message and exit
  -l LOVD, --lovd LOVD  Input LOVD file
  --region REGION       Input region file
  -p PADDING, --padding PADDING
                        Padding around regions or variants in bases. Defaults
                        to 100
  -x XML, --xml XML     Output Miracle XML file
  -t TSV, --tsv TSV     Output TSV file
  -b BAM, --bam BAM     Path to output BAM file containing primers
  -s SAMPLE, --sample SAMPLE
                        Same ID for regions
  --product_size PRODUCT_SIZE
                        Size range of desired product. Defaults to 200-450
                        This will be taken as a minimum product size in the
                        case of regions
  --min-margin MIN_MARGIN
                        Minimum distance from region or variant. Default = 10
  --strict              Enable strict mode. Primers with products larger than
                        max product size will NOT be returned
  --n_raw_primers N_RAW_PRIMERS
                        Legacy option. Will be ignored
  --m13                 Output primers with m13 tails
  --m13-forward M13_FORWARD
                        Sequence of forward m13 tail
  --m13-reverse M13_REVERSE
                        Sequence of reverse m13 tail
  -f FIELD, --field FIELD
                        Name of field in DBSNP file for allele frequency
  -af ALLELE_FREQ, --allele-freq ALLELE_FREQ
                        Max accepted allele freq
  -fq1 FQ1              Path to forward fastq file for primer output. Set if
                        you want to export your primers in fastq format
                        (qualities will be sanger-encoded 40)
  -fq2 FQ2              Path to reverse fastq file for primer output. Set if
                        you want to export your primers in fastq format
                        (qualities will be sanger-encoded 40)
  -R REFERENCE, --reference REFERENCE
                        Path to reference fasta file
  --dbsnp DBSNP         Path to DBSNP vcf
  --primer3 PRIMER3     Path to primer3_core exe
  --bwa BWA             Path to BWA exe
  --samtools SAMTOOLS   Path to samtools exe
  --ignore-errors       Ignore errors
  --opt-primer-length OPT_PRIMER_LENGTH
                        Optimum primer length (default = 25)
  --opt-gc-perc OPT_GC_PERC
                        Optimum primer GC percentage (default = 50)
  --min-melting-temperature MIN_MELTING_TEMPERATURE
                        Minimum primer melting temperature (default = 58)
  --max-melting-temperature MAX_MELTING_TEMPERATURE
                        Maximum primer melting temperature (default = 62)

```

### Usage

`primerdesign` requires either the `-l` or `--region` arguments for input.

* `-l` : LOVD import file format
* `--region`: BED track
 
 When `--region` is used, a sample name must be given with `-s <sample>`. 

It will output in either of the following two formats with the `-x` and `-t` arguments:

* `-x`: Miracle XML
* `-t`: BED-like tab-delimited format

The following arguments are mandatory:

* `-R`: path to reference fasta. This fasta *must* have a `.fai` index. Generate this with `samtools faidx <fasta>`
* `--dbsnp`: Path to DBSNP VCF file
* `--primer3`: Path to primer 3 executable.

Recommended arguments:

* `--samtools`: Path to samtools. If not given, will simply assume `samtools` is on the PATH.
* `--bwa`: Path to bwa. If not given, will simply assume `bwa` is on the PATH. 


Known issues
------------
1. In some cases, no suitable primer pair can be found. If `primerdesign` fails with a `NoPrimersException` you can try to increase the padding.
 
 