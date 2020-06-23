GDC Variant Filtration Tool
---

This repository contains the source code used in the VCF variant filtration
workflows within the GDC. A single CLI is generated with multiple subcommands. 

## Requirements 

* Python >= 3.6
* pysam
* defopt
 
## Subcommands 

### `add-oxog-filters`

Adds 'oxog' filter tag to VCFs.

```
usage: gdc-filtration-tools add-oxog-filters [-h]
                                             input_vcf input_dtoxog output_vcf

Adds 'oxog' filter tag to VCFs.

positional arguments:
  input_vcf     The full input VCF file to filter.
  input_dtoxog  The dtoxog VCF from dtoxog-maf-to-vcf used to annotate the full input VCF.
  output_vcf    The output filtered VCF file to create.     BGzip and tabix-index created if ends with '.gz'.

optional arguments:
  -h, --help    show this help message and exit
```

### `create-dtoxog-maf`

Takes a SNP-only VCF file and converts it to the dToxoG MAF format
which includes the OXOQ value.

```
usage: gdc-filtration-tools create-dtoxog-maf [-h]
                                              input_vcf output_file reference
                                              oxog_file oxoq_score

Takes a SNP-only VCF file and converts it to the dToxoG MAF format
which includes the OXOQ value.

positional arguments:
  input_vcf    The input SNP-only VCF file to convert to dToxoG MAF.
  output_file  The output MAF file to create.
  reference    Faidx indexed reference fasta file.
  oxog_file    Metrics file output from GATK OxoGMetrics tool.
  oxoq_score   The oxoQ score.

optional arguments:
  -h, --help   show this help message and exit
```

### `create-oxog-intervals`

Takes a SNP-only VCF file and creates an interval list for
use by the Broad oxog metrics tool.

```
usage: gdc-filtration-tools create-oxog-intervals [-h] input_vcf output_file

Takes a SNP-only VCF file and creates an interval list for
use by the Broad oxog metrics tool.

positional arguments:
  input_vcf    The input SNP-only VCF file to extract intervals from.
  output_file  The output interval list to create.

optional arguments:
  -h, --help   show this help message and exit
```

### `dtoxog-maf-to-vcf`

Transforms dToxoG MAF to minimal VCF of only dtoxo failures.

```
usage: gdc-filtration-tools dtoxog-maf-to-vcf [-h]
                                              input_maf reference_fa
                                              output_vcf

Transforms dToxoG MAF to minimal VCF of only dtoxo failures.

positional arguments:
  input_maf     The annotated dtoxog MAF output file.
  reference_fa  Reference fasta used to make seqdict header.
  output_vcf    The output minimal VCF with only failed dtoxog records     BGzip and tabix-index created if ends with '.gz'.

optional arguments:
  -h, --help    show this help message and exit
```

### `extract-oxoq-from-sqlite`

Extract the OXOQ score for a particular context from the GDC
harmonization metrics SQLite file. The score is printed to
stdout.

```
usage: gdc-filtration-tools extract-oxoq-from-sqlite [-h] [-c CONTEXT]
                                                     [-t TABLE]
                                                     [-i INPUT_STATE]
                                                     db_file

Extract the OXOQ score for a particular context from the GDC
harmonization metrics SQLite file. The score is printed to
stdout.

positional arguments:
  db_file               Path to the SQLite db file.

optional arguments:
  -h, --help            show this help message and exit
  -c CONTEXT, --context CONTEXT
                        The nucleotide context of interest.
                        (default: CCG)
  -t TABLE, --table TABLE
                        The SQLite table name.
                        (default: picard_CollectOxoGMetrics)
  -i INPUT_STATE, --input-state INPUT_STATE
                        The input state to select for in the input_state column.
                        (default: markduplicates_readgroups)
```

### `filter-contigs`

Filter out VCF records on chromosomes that are not present
in the contig lines of the VCF header.

```
usage: gdc-filtration-tools filter-contigs [-h] input_vcf output_vcf

Filter out VCF records on chromosomes that are not present
in the contig lines of the VCF header.

positional arguments:
  input_vcf   The input VCF file to filter.
  output_vcf  The output filtered VCF file to create.     BGzip and tabix-index created if ends with '.gz'.

optional arguments:
  -h, --help  show this help message and exit
```

### `filter-nonstandard-variants`

Remove non-ACTG loci from a SNP-ONLY VCF. No validation that
the VCF is SNP-only is done.

```
usage: gdc-filtration-tools filter-nonstandard-variants [-h]
                                                        input_vcf output_vcf

Remove non-ACTG loci from a SNP-ONLY VCF. No validation that
the VCF is SNP-only is done.

positional arguments:
  input_vcf   The input SNP-only VCF file to filter.
  output_vcf  The output filtered VCF file to create.     BGzip and tabix-index created if ends with '.gz'.

optional arguments:
  -h, --help  show this help message and exit
```

### `filter-somatic-score`

Filters SomaticSniper VCF files based on the Somatic Score.

```
usage: gdc-filtration-tools filter-somatic-score [-h] [-t TUMOR_SAMPLE_NAME]
                                                 [-d DROP_SOMATIC_SCORE]
                                                 [-m MIN_SOMATIC_SCORE]
                                                 input_vcf output_vcf

Filters SomaticSniper VCF files based on the Somatic Score.

positional arguments:
  input_vcf             The input VCF file to filter.
  output_vcf            The output filtered VCF file to create.     BGzip and tabix-index created if ends with '.gz'.

optional arguments:
  -h, --help            show this help message and exit
  -t TUMOR_SAMPLE_NAME, --tumor-sample-name TUMOR_SAMPLE_NAME
                        The name of the tumor sample in the VCF.
                        (default: TUMOR)
  -d DROP_SOMATIC_SCORE, --drop-somatic-score DROP_SOMATIC_SCORE
                        If the somatic score is < this, remove it.
                        (default: 25)
  -m MIN_SOMATIC_SCORE, --min-somatic-score MIN_SOMATIC_SCORE
                        If the somatic score is > drop_somatic_score                               and < this value, add ssc filter tag.
                        (default: 40)
```

### `format-gdc-vcf`

Adds VCF header metadata specific to the GDC.

```
usage: gdc-filtration-tools format-gdc-vcf [-h] [-r REFERENCE_NAME]
                                           input_vcf output_vcf
                                           patient_barcode case_id
                                           tumor_barcode tumor_aliquot_uuid
                                           tumor_bam_uuid normal_barcode
                                           normal_aliquot_uuid normal_bam_uuid

Adds VCF header metadata specific to the GDC.

positional arguments:
  input_vcf             The input VCF file to format.
  output_vcf            The output formatted VCF file to create. BGzip and tabix-index created if ends with '.gz'.
  patient_barcode       The case submitter id.
  case_id               The case uuid.
  tumor_barcode         The tumor aliquot submitter id.
  tumor_aliquot_uuid    The tumor aliquot uuid.
  tumor_bam_uuid        The tumor bam uuid.
  normal_barcode        The normal aliquot submitter id.
  normal_aliquot_uuid   The normal aliquot uuid.
  normal_bam_uuid       The normal bam uuid.

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE_NAME, --reference-name REFERENCE_NAME
                        Reference name to use in header.
                        (default: GRCh38.d1.vd1.fa)
```

### `format-pindel-vcf`

Formats Pindel VCFs to work better with GDC downstream workflows.

```
usage: gdc-filtration-tools format-pindel-vcf [-h] input_vcf output_vcf

Formats Pindel VCFs to work better with GDC downstream workflows.

positional arguments:
  input_vcf   The input VCF file to filter.
  output_vcf  The output filtered VCF file to create. BGzip and tabix-index created if ends with '.gz'.

optional arguments:
  -h, --help  show this help message and exit
```

### `position-filter-dkfz`

Removes VCF records where the POS-2 is less than 0 which
will cause an Exception to be thrown in DKFZBiasFilter. We
assume that the input VCF only contains SNPs, but no assertions
are made to validate this.

```
usage: gdc-filtration-tools position-filter-dkfz [-h] input_vcf output_vcf

Removes VCF records where the POS-2 is less than 0 which
will cause an Exception to be thrown in DKFZBiasFilter. We
assume that the input VCF only contains SNPs, but no assertions
are made to validate this.

positional arguments:
  input_vcf   The input VCF file to filter.
  output_vcf  The output filtered VCF file to create. BGzip and tabix-index created if ends with '.gz'.

optional arguments:
  -h, --help  show this help message and exit
```

## Docker Tools

**variant-filtration-tool** <br />
Repository: https://quay.io/ncigdc/variant-filtration-tool <br />
Tag: 0.1 <br />
