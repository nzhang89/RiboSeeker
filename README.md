# RiboSeeker: An End-to-End Package for Ribosome Profiling Data Analysis

Ribosome profiling is a technology for determining translation activity by sequencing ribosome protected mRNA fragments. We developed RiboSeeker, an end-to-end package for ribosome profiling data analysis. RiboSeeker has two components: [Ribomake](https://github.com/nzhang89/Ribomake) and [RiboSeeker](https://github.com/nzhang89/RiboSeeker):

* Ribomake: a Snakemake workflow to process ribosome profiling data from raw reads;
* RiboSeeker: an R package for further data analysis from aligned reads;

Please refer to [Ribomake](https://github.com/nzhang89/Ribomake) page for more details about the installation process and usage. This repository focuses on the R package to process aligned Ribo- and RNA-seq reads.

### Installation

You can install the RiboSeeker R package from [GitHub](https://github.com/nzhang89/RiboSeeker) with

```r
if (!requireNamespace('devtools', quietly = TRUE)) {
  install.packages('devtools')
}
devtools::install_github('nzhang89/RiboSeeker')
```

### Usage

We provide a simple [tutorial]() to demonstrate the basic usage of the RiboSeeker R package.
