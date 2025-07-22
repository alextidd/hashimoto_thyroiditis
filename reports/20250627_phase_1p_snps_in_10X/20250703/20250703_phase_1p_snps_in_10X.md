---
title: "Phasing 1p SNPs in PD63118 and genotyping in 10X"
author: "Alexandra Tidd"
date: "03 July, 2025"
output:
  html_document:
    fig_width: 8
    keep_md: true
    toc: true
    toc_float: true
    toc_collapsed: true
    toc_depth: 4
    theme: lumen
---

The longest 1p LOH event is observed in cells `plate3_wellA2`, affecting the
entire 1p arm, and is therefore expected to be homozygous for all SNPs. We can
therefore phase all SNPs in PD63118 chr1p. Then, I will look for genic SNPs, as 
these will be most likely to be picked up in the snRNAseq.



We get the cell IDs with 1p LOH.



First, we get all SNPs on 1p in the LOH cells.





Next, we phase the SNPs on 1p in `plate3_wellA2`.



Now, we visualise the phased SNPs for `plate3_wellA2` in a BAF plot.

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_baf_plate3_wellA2_snps-1.png" style="display: block; margin: auto;" />

## Phased BAF plots

Now we plot this phasing in all other LOH cell.

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-1.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-2.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-3.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-4.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-5.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-6.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-7.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-8.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-9.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-10.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-11.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-12.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-13.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-14.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-15.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-16.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-17.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-18.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-19.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-20.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-21.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-22.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-23.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-24.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-25.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-26.png" style="display: block; margin: auto;" /><img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/plot_bafs-27.png" style="display: block; margin: auto;" />

## Breakpoints

Next, we attempt to identify the breakpoints in each case.


```
## [1] "plate10_wellA10_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-1.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellA8_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-2.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellA9_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-3.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellC3_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-4.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellD2_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-5.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellD3_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-6.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellD6_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-7.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellD7_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-8.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellF12_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-9.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellG10_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-10.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellG2_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-11.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellG4_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-12.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellG7_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-13.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellG9_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-14.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellH4_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-15.png" style="display: block; margin: auto;" />

```
## [1] "plate10_wellH5_dna_run50382"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-16.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellA2_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-17.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellB3_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-18.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellC8_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-19.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellD4_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-20.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellD7_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-21.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellE10_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-22.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellE11_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-23.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellF6_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-24.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellF7_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-25.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellH5_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-26.png" style="display: block; margin: auto;" />

```
## [1] "plate3_wellH6_dna_run49882"
```

<img src="/nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/reports/20250627_phase_1p_snps_in_10X_files/figure-html/id_breakpoints-27.png" style="display: block; margin: auto;" />
