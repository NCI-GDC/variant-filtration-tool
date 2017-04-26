GDC Variant Filtration Tool
---

This tool only runs the somatic score (somatic sniper only) filtration

## Inputs

* VCF file you want to filter

## Outputs

* Filtered VCF file

## Usage

```
$ python main.py somaticscore -h
usage: main.py somaticscore [-h] --input_vcf INPUT_VCF --output_vcf OUTPUT_VCF
                            [--tumor_sample_name TUMOR_SAMPLE_NAME]
                            [--drop_somatic_score DROP_SOMATIC_SCORE]
                            [--min_somatic_score MIN_SOMATIC_SCORE]

optional arguments:
  -h, --help            show this help message and exit
  --input_vcf INPUT_VCF
                        The VCF you want to filter
  --output_vcf OUTPUT_VCF
                        The output filtered VCF file. Will gzip if ends in
                        ".gz"
  --tumor_sample_name TUMOR_SAMPLE_NAME
                        The name of the TUMOR sample in the VCF
  --drop_somatic_score DROP_SOMATIC_SCORE
                        If the somatic score is < this, remove it [25]
  --min_somatic_score MIN_SOMATIC_SCORE
                        If the somatic score is > drop_somatic_score and <
                        this value, add ssc filter tag [40]
```

## Docker Tools

**variant-filtration-tool** <br />
Repository: https://quay.io/ncigdc/variant-filtration-tool <br />
Tag: 0.1 <br />
Description: Contains fpfilter, bamreadcount, samtools, this script, and requirements
