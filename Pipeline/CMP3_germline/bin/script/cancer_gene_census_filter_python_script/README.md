# cancer_gene_census_filter.py
A Python script using a file containing cancer genes to filter a VCF file. 
The file cancer_gene_census.csv, available at the link: https://cancer.sanger.ac.uk/cosmic/download, can be used to select from the vcf file genes which contain mutations that have been causally implicated in cancer.

## Usage
From the directory containing the cancer_gene_census_filter.py script launch it with all the required parameters.

Example:
```bash
python3 cancer_gene_census_filter.py -v vcf_to_be_fileted.vcf -g file_of_gene.csv

```

## Needed parameters
* **-v**: the path to the VCF file to be filtered.
* **-g**: the path to the file containig genes to be selected from the VCF.
* **-s**: the separator of columns in the genes file (DEFAULT=",")
* **-c**: the column selected from the genes file (DEFAULT=1)
* **-e**: this option exclude from the vcf the genes included.
* **-a**: annotator chooser, s for snpeff and v for vep (DEFAULT=s). 


License 
[MIT](https://choosealicense.com/licenses/mit/)
