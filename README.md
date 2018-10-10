PRINIA
=======

Prinia is a python package for designing primers from genomic regions
or LOVD import format files.

Requirements
-------------

## Python version

3.5+

### Python packages

* future >= 0.15.2
* biopython >= 1.66
* lxml >= 3.4.4
* pyfaidx >= 0.4.4
* pyvcf >= 0.6.7
* jsonschema >= 2.6.0
* Jinja2 >= 2.9.5

## Other dependencies

You must have a recent version of [primer3](http://primer3.sourceforge.net/), [samtools](http://www.htslib.org/download/)
and [bwa](https://github.com/lh3/bwa)


Installation
-------------

It is recommended you use a [virtual environment](https://virtualenv.readthedocs.org),
such that the risk of dependency hell is minimized.
  
After creating your virtualenv, clone this repository. 
Install all the required python packages with:
```
pip install -r requirements.txt
```


Then install PRINIA with:

```
python setup.py install
```


Primer design
-------------

A `primerdesign` tool will be be added to your path upon installation
of prinia. This tool will generate your primers from BED records
(regions) or LOVD import file format files.
 
### Help


```
usage: primerdesign [-h] (-l LOVD | --region REGION) [-p PADDING]
                    (-x XML | -t TSV) -b BAM [-s SAMPLE]
                    [--min-margin MIN_MARGIN] [--strict] [--m13]
                    [--m13-forward M13_FORWARD] [--m13-reverse M13_REVERSE]
                    [-f FIELD] [-af ALLELE_FREQ] [-fq1 FQ1] [-fq2 FQ2] -R
                    REFERENCE [--dbsnp DBSNP] --primer3 PRIMER3 [--bwa BWA]
                    [--samtools SAMTOOLS] [--settings-json SETTINGS_JSON]
                    [--ignore-errors]

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
  --min-margin MIN_MARGIN
                        Minimum distance from region or variant. Default = 10
  --strict              Enable strict mode. Primers with products larger than
                        max product size will NOT be returned
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
  --settings-json SETTINGS_JSON
                        Optional path to primer3 settings json file.
  --ignore-errors       Ignore errors
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


### Primer3 settings

Some settings that are passed on to primer3 can be provided with a json
file. Said json file must conform to the json schema provided
[here](prinia/static/primer3_settings_schema.json). A validation error
is thrown if the provided json does not conform to the schema.

Parameters not given in the settings file are supplied default values.
In case no settings file is given at all, only the default values are
used.

The schema defines all default values. For convenience, they are listed
here as well:

| parameter | explanation | default |
| --------- | ----------- | ------- |
| `primer_min_gc` | Minimum GC-percentage of primers | 20 |
| `primer_internal_min_gc` | Equivalent parameter of primer_min_gc for the internal oligo. | 20 |
| `primer_opt_gc_percent` | Optimum GC-percentage of primers | 50 |
| `primer_max_gc` | Maximum GC-percentage of primers. | 80 |
| `primer_internal_max_gc` | Equivalent parameter of primer_max_gc for the internal oligo. | 80 |
| `primer_wt_gc_percent_lt` | Penalty weight for primers with GC-percentage lower than primer_opt_gc_percent | 0 |
| `primer_internal_wt_gc_percent_lt` | Equivalent parameter of primer_wt_gc_percent_lt for the interal oligo. | 0 |
| `primer_wt_gc_percent_gt` | Penalty weight for primers with GC-percentage higher than primer_opt_gc_percent | 0 |
| `primer_internal_wt_gc_percent_gt` | Equivalent parameter of primer_wt_gc_percent_gt for the internal oligo. | 0 |
| `primer_gc_clamp` | Require the specified number of consecutive Gs and Cs at the 3' end of both the left and right primer. | 0 |
| `primer_max_end_gc` | The maximum number of Gs or Cs allowed in the last five 3' bases of a left or right primer. | 5 |
| `primer_opt_size` | Optimum size of primers in bases. | 25 |
| `primer_min_size` | Minimum size of primers in bases.  | 20 |
| `primer_max_size` | Maximum size of primers in bases. | 30 |
| `primer_max_ns_accepted` | Maximum number of unknown bases (N) accepted in a primer | 0 |
| `primer_product_size_range` | Accepted size range of the product size. | 200-450 |
| `primer_product_opt_size` | Optimum size of product. Will default to the midpoint in primer_product_size_range unless specified | 325 |
| `primer_pair_wt_product_size_gt` | Penalty weight for products larger than primer_product_opt_size | 0.1 |
| `primer_pair_wt_product_size_lt` | Penalty weight for products smaller than primer_product_opt_size | 0.1 |
| `primer_min_tm` | Minimum melting temperature of primer degrees Celsius | 58 |
| `primer_max_tm` | Maximum melting temperature in degrees Celsius | 62 |
| `primer_num_return` | Number of returned primers. Increase to increase search space | 200 |


Some example configuration files can be found [here](tests/data/valid_settings.json), [here](tests/data/valid_partial_settings.json) and [here](tests/data/valid_partial_settings2.json).

Known issues
------------
1. In some cases, no suitable primer pair can be found. If `primerdesign` fails with a `NoPrimersException` you can try to increase the padding.
 
 