## Setup

- Set-up: [R and RStudio](https://uclouvain-cbio.github.io/WSBIM1207/sec-startr.html)
- [Bioconductor](https://uclouvain-cbio.github.io/WSBIM1207/sec-bioinfo.html)
- [Objects](https://uclouvain-cbio.github.io/WSBIM1322/sec-obj.html)

## Introduction

- [How does mass spectrometry work?](https://lgatto.github.io/bioc-ms-prot/lab.html#2_how_does_mass_spectrometry_work)
- [Accessing raw data](https://lgatto.github.io/bioc-ms-prot/lab.html#3_accessing_data)

## Raw data

## Identification data

## Quantitative data

Pre-req: `BiocManager::install("statOmics/MSqRob")`

- https://lgatto.github.io/bioc-ms-prot/lab.html#6_quantitative_proteomics
- Processing: https://lgatto.github.io/bioc-ms-prot/csama2019-lab.html#quantitative-data
- Missing values: https://lgatto.github.io/bioc-ms-prot/lab.html#missing_values
- [Stats intro](https://lgatto.github.io/bioc-ms-prot/lab.html#8_statistical_analysis)


```{r}
## using limma
design <- model.matrix(~ cptac_rob$condition)
fit <- lmFit(exprs(cptac_rob), design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "cptac_rob$conditionB", number = Inf)

library("tidyverse")
fData(cptac_rob) <-
    full_join(rownames_to_column(fData(cptac_rob)),
              rownames_to_column(res)) %>%
    column_to_rownames()

library("ggrepel")
ggplot(fData(cptac_rob),
       aes(x = logFC,
           y = -log10(adj.P.Val))) +
    geom_point()

ggplot(fData(cptac_rob),
       aes(x = logFC,
           y = -log10(adj.P.Val),
           label = sub("\\[OS.+\\]", "", Proteins))) +
    geom_point() +
    geom_text_repel(
        data = subset(fData(cptac_rob), adj.P.Val < 0.05),
        nudge_x = 0.05,
        nudge_y = -0.05,
        segment.size = 0.5,
    )

```

- [Differential expression with `MSqRob`](https://lgatto.github.io/bioc-ms-prot/csama2019-lab.html#differential-expression-analysis)

## R for mass spectrometry

The aim of the [R for mass
spectrometry](https://www.rformassspectrometry.org/)
RforMassSpectrometry initiative is to provide efficient, thoroughly
documented, tested and flexible R software for the analysis and
interpretation of high throughput mass spectrometry assays, including
proteomics and metabolomics experiments. The project formalises the
longtime collaborative development efforts of its core members under
the RforMassSpectrometry organisation to facilitate dissemination and
accessibility of their work.

Two packages of interest:

- [Specra](https://rformassspectrometry.github.io/Spectra/)
- [Features](https://rformassspectrometry.github.io/Features/)


## Use cases

- (Differential expression)
- Protein-protein interaction (BioID2)

