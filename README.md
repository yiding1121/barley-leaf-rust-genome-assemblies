# Barley leaf rust genome assemblies

This repository contains genome annotations and analysis code for haplotype-resolved assemblies of the barley leaf rust fungus *Puccinia hordei*. The study compares isolates Ph518 and Ph560 and focuses on chromosome-wide comparisons, genome structure, population variation, copy-number variation (CNV), effector expression, and the example gene *Cyp51*.

## Analysis workflow

1. Assemble and annotate both haplotypes of Ph518 and Ph560.
2. Assess chromosome-scale assembly structure using Hi-C and whole-genome alignments.
3. Compare chromosome 9 gene, effector, transposable-element, and Hi-C features.
4. Analyze population structure using SNPs and chromosome 9 k-mers.
5. Identify CNVs and genes affected by deletions or duplications.
6. Examine population-genetic statistics, effector expression, and the *Cyp51* locus.

## Files

| File | Analysis |
| --- | --- |
| [`Ph518A_braker.EVM.gff3.gz`](Ph518A_braker.EVM.gff3.gz), [`Ph518B_braker.EVM.gff3.gz`](Ph518B_braker.EVM.gff3.gz) | Ph518 gene annotations |
| [`Ph560A_braker.EVM.gff3.gz`](Ph560A_braker.EVM.gff3.gz), [`Ph560B_braker.EVM.gff3.gz`](Ph560B_braker.EVM.gff3.gz) | Ph560 gene annotations |
| [`Chr9_composition.R`](Chr9_composition.R) | Chromosome 9 genes, effectors, and structural breakpoints |
| [`TE_and_HiC_analysis.R`](TE_and_HiC_analysis.R) | Gene density, TE density, GC content, and Hi-C features |
| [`Chr9_composition_and_statistics.R`](Chr9_composition_and_statistics.R) | Chromosome 9 composition and statistical comparisons |
| [`Population_GOenrichment_and_Pi_TajimasD_Fst.R`](Population_GOenrichment_and_Pi_TajimasD_Fst.R) | Population structure, GO enrichment, nucleotide diversity, Tajima's D, and FST |
| [`k-mer_analysis.ipynb`](k-mer_analysis.ipynb) | Chromosome 9 k-mer clustering |
| [`CNV_analysis_full.R`](CNV_analysis_full.R) | CNV detection and gene-overlap analysis |
| [`effector_expression_chr9.R`](effector_expression_chr9.R) | Effector expression and TE-proximity analysis |

> The scripts were developed interactively and contain paths to data files that are not stored in this repository. Update the input paths before running them.

## Citation

Please cite the associated barley leaf rust genome study when using these assemblies, annotations, or analysis code:

Yu, X., Chhetri, M., Hasan, M., Roberts, T.H., Zhang, P., Park, R.F. and Ding, Y. (2026), Comparative genomics reveals lineage-associated structural variation and diversification in a barley fungal pathogen. New Phytol. https://doi.org/10.1111/nph.71571


