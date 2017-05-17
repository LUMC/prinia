PRINIA
=======

Prinia is a python package for designing primers from genomic regions or LOVD import format files. 

Requirements
-------------

You must have a recent version of [primer3](http://primer3.sourceforge.net/) and [blat](http://sourceforge.net/projects/blat/files/Blat%20Full%20Version/)

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
                    (-x XML | -t TSV) [-s SAMPLE]
                    [--product_size PRODUCT_SIZE]
                    [--n_raw_primers N_RAW_PRIMERS] [--m13]
                    [--m13-forward M13_FORWARD] [--m13-reverse M13_REVERSE]
                    [-f FIELD] [-af ALLELE_FREQ] -R REFERENCE --dbsnp DBSNP
                    --primer3 PRIMER3 --blat BLAT

optional arguments:
  -h, --help            show this help message and exit
  -l LOVD, --lovd LOVD  Input LOVD file
  --region REGION       Input region file
  -p PADDING, --padding PADDING
                        Padding around regions or variants in bases. Defaults
                        to 100
  -x XML, --xml XML     Output Miracle XML file
  -t TSV, --tsv TSV     Output TSV file
  -s SAMPLE, --sample SAMPLE
                        Same ID for regions
  --product_size PRODUCT_SIZE
                        Size range of desired product. Defaults to 200-450
                        This will be taken as a minimum product size in the
                        case of regions
  --n_raw_primers N_RAW_PRIMERS
                        Amount of raw primers from primer3 output that will be
                        considered. By default, only the top 4 primers (in
                        each direction) will be considered. Increasing this
                        value does not mean more primers will be returned; it
                        simply increases the search space.
  --m13                 Output primers with m13 tails
  --m13-forward M13_FORWARD
                        Sequence of forward m13 tail
  --m13-reverse M13_REVERSE
                        Sequence of reverse m13 tail
  -f FIELD, --field FIELD
                        Name of field in DBSNP file for allele frequency
  -af ALLELE_FREQ, --allele-freq ALLELE_FREQ
                        Max accepted allele freq
  -R REFERENCE, --reference REFERENCE
                        Path to reference fasta file
  --dbsnp DBSNP         Path to DBSNP vcf
  --primer3 PRIMER3     Path to primer3 exe
  --blat BLAT           Path to blat exe

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
* `--primer3`: Path to primer 3 executable. On the SHARK cluster, this is `/usr/bin/primer3_core`
* `--blat`: Path to blat executable. On the SHARK cluster, this is `/usr/local/bin/blat` 


Known issues
------------
1. In some cases, no suitable primer pair can be found. If `primerdesign` fails with a `NoPrimersException` you can try to increase the padding.
 
 