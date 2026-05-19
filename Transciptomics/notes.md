
### Gene annotations.
```bash
mkdir gtf
cd gtf
wget ftp://ftp.ensembl.org/pub/release-112/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.112.gtf.gz
gunzip Caenorhabditis_elegans.WBcel235.112.gtf.gz

cd ..
Rscript gene_annotations.R
```

### RNA-seq two-group comparisons – differential expression at each age.
```bash
Rscript two_group_comparisons.R
Rscript two_group_volcano_plots.R
```

### Time course comparison – MDS plot of all samples.
```bash
Rscript time_course_comparison.R
```

### Enrichment analysis and plots.
```bash
Rscript enrichment_analysis.R
Rscript enrichment_plots.R
````
